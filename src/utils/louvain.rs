use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
use ordered_float::NotNan;
use std::{
    cmp::Reverse,
    env::args,
    fmt,
    fs::File,
    io::{BufRead, BufReader},
    str::FromStr,
};

#[derive(Debug, Default)]
pub struct Node {
    pub id: u32,             // community ID
    pub weight: f32,         //weight sum between children nodes
    pub nodes: HashSet<u32>, //the children nodes id in this parent node
}

impl Node {
    fn new(id: u32, weight: f32, nodes: impl IntoIterator<Item = u32>) -> Self {
        Node {
            id,
            weight,
            nodes: HashSet::from_iter(nodes),
        }
    }
}

#[derive(Debug, Default)]
pub struct Louvain {
    data: HashMap<u32, HashMap<u32, f32>>,   //raw data and weight
    communities: HashMap<u32, HashSet<u32>>, // community ID, which contains the set of Node
    node: HashMap<u32, Node>,                // Node and its instance
}

impl fmt::Display for Louvain {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (id, nodes) in self.communities.iter().filter(|(_, v)| !v.is_empty()) {
            write!(
                f,
                "id: {} count: {} eles: ",
                id,
                nodes
                    .iter()
                    .map(|x| self.node[x].nodes.len())
                    .sum::<usize>()
            )?;
            for node in nodes {
                for n in &self.node[node].nodes {
                    write!(f, "{} ", n)?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl Louvain {
    pub fn new(data: HashMap<u32, HashMap<u32, f32>>) -> Self {
        let mut lv = Louvain {
            data,
            ..Default::default()
        };
        for vid in lv.data.keys().copied() {
            lv.communities.insert(vid, HashSet::from_iter([vid]));
            lv.node.insert(vid, Node::new(vid, 0., [vid]));
        }
        lv
    }

    fn first_stage(&mut self) -> bool {
        let mut mod_inc = false;
        let mut node_ids: HashMap<u32, f32> = HashMap::default();
        let mut visit_ids: Vec<u32> = self.data.keys().copied().collect();

        visit_ids.sort(); //here can make the speed increase twice as fast
        loop {
            //safe unwrap
            let mut can_stop = true;
            for v_id in &visit_ids {
                let v_nid = self.node[v_id].id;
                node_ids.clear();

                for w_nid in self.data[v_id].keys().map(|x| self.node[x].id) {
                    if node_ids.contains_key(&w_nid) {
                        continue;
                    }
                    let communities = &self.communities[&w_nid];
                    node_ids.insert(
                        w_nid,
                        self.data[v_id]
                            .iter()
                            .filter_map(|(k, v)| communities.get(k).map(|_| v))
                            .sum(),
                    );
                }

                if let Some((id, max_weight)) = node_ids
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap().then(b.0.cmp(a.0)))
                {
                    if *max_weight > 0.0 && *id != v_nid {
                        self.node.get_mut(v_id).unwrap().id = *id;
                        self.communities.get_mut(id).unwrap().insert(*v_id);
                        self.communities.get_mut(&v_nid).unwrap().remove(v_id);
                        can_stop = false;
                        mod_inc = true;
                    }
                }
            }
            if can_stop {
                break;
            }
        }
        mod_inc
    }

    fn second_stage(mut self) -> Self {
        let mut node: HashMap<u32, Node> = HashMap::default();
        let mut communities: HashMap<u32, HashSet<u32>> = HashMap::default();
        let mut decluster_ids = Vec::new();
        for (id, nodes) in self.communities.iter().filter(|(_, v)| !v.is_empty()) {
            let mut new_node = Node::new(*id, 0., []);
            for nid in nodes {
                let vertex = &self.node[nid];
                new_node.nodes.extend(&vertex.nodes);
                new_node.weight += vertex.weight;
                if let Some(nid_v) = self.data.get(nid) {
                    for (_, v) in nid_v.iter().filter(|(k, _)| nodes.contains(k)) {
                        new_node.weight += v / 2.0;
                    }
                }
            }

            if new_node.weight < 0. {
                decluster_ids.push(*id);
            } else {
                communities.insert(*id, HashSet::from_iter([*id]));
                node.insert(*id, new_node);
            }
        }

        //decluster communities with weight < 0.
        for id in decluster_ids {
            let nodes = self.communities.remove(&id).unwrap(); //safely
            for nid in nodes {
                // if communities or node already contain nid
                let mut new_nid = nid;
                while communities.contains_key(&new_nid) || node.contains_key(&new_nid){
                    new_nid += 1;
                }

                communities.insert(new_nid, HashSet::from_iter([new_nid]));
                node.insert(
                    new_nid,
                    Node::new(
                        new_nid,
                        self.node[&nid].weight,
                        self.node[&nid].nodes.iter().copied(),
                    ),
                );
                self.communities.insert(new_nid, HashSet::from_iter([nid]));
            }
        }

        let mut data: HashMap<u32, HashMap<u32, f32>> = HashMap::default();
        for (nid1, nodes1) in self.communities.iter().filter(|(_, v)| !v.is_empty()) {
            for (nid2, nodes2) in self
                .communities
                .iter()
                .filter(|(k, v)| *k > nid1 && !v.is_empty())
            {
                let mut edge_weight = 0.0;
                for vid in nodes1 {
                    if let Some(vid_v) = self.data.get(vid) {
                        for (_, v) in vid_v.iter().filter(|(k, _)| nodes2.contains(k)) {
                            edge_weight += v;
                        }
                    }
                }

                if edge_weight != 0. {
                    insert_data(&mut data, *nid1, *nid2, edge_weight);
                    insert_data(&mut data, *nid2, *nid1, edge_weight);
                }
            }
        }

        Self {
            data,
            communities,
            node,
        }
    }

    fn get_communities(self) -> (HashMap<u32, HashMap<u32, f32>>, Vec<Node>) {
        let mut communities = Vec::new();
        for (&id, nodes) in self.communities.iter().filter(|(_, v)| !v.is_empty()) {
            let mut weight = 0.;
            let mut new_nodes = HashSet::default();
            for vid in nodes {
                let v = &self.node[vid];
                new_nodes.extend(&v.nodes);
                weight += v.weight; //weight of v
                if let Some(ks) = self.data.get(vid) {
                    //weight between v
                    for (_, v) in ks.iter().filter(|(k, _)| nodes.contains(k)) {
                        weight += v / 2.0;
                    }
                }
            }
            communities.push(Node {
                id,
                weight,
                nodes: new_nodes,
            });
        }

        //init a new data
        let mut data: HashMap<u32, HashMap<u32, f32>> = HashMap::default();
        for c1 in &communities {
            for c2 in communities.iter().filter(|x| x.id > c1.id) {
                let mut weight = 0.;
                for n1 in &self.communities[&c1.id] {
                    for n2 in &self.communities[&c2.id] {
                        weight += self
                            .data
                            .get(n1)
                            .map_or_else(|| 0., |v| *v.get(n2).unwrap_or(&0.));
                    }
                }
                if weight != 0. {
                    assert!(
                        weight < 0.,
                        "the weight of two conflicting community is not less than 0"
                    );
                    insert_data(&mut data, c1.id, c2.id, weight);
                    insert_data(&mut data, c2.id, c1.id, weight);
                }
            }
        }

        (data, communities)
    }

    pub fn execute(mut self) -> (HashMap<u32, HashMap<u32, f32>>, Vec<Node>) {
        loop {
            let mod_inc = self.first_stage();
            if mod_inc {
                self = self.second_stage();
            } else {
                return self.get_communities();
            }
        }
    }
}

pub fn new_data() -> HashMap<u32, HashMap<u32, f32>> {
    HashMap::default()
}

#[allow(dead_code)]
fn out_data(data: &HashMap<u32, HashMap<u32, f32>>) {
    for (k, v) in data.iter() {
        for (v, w) in v.iter() {
            println!("G: {k}\t{v}\t{w}");
        }
    }
}

#[inline(always)]
pub fn insert_data(data: &mut HashMap<u32, HashMap<u32, f32>>, k1: u32, k2: u32, v: f32) {
    data.entry(k1)
        .and_modify(|x| {
            x.entry(k2).and_modify(|x| *x += v).or_insert(v);
        })
        .or_insert_with(|| HashMap::from_iter([(k2, v)]));
}

#[inline(always)]
pub fn assign_data(data: &mut HashMap<u32, HashMap<u32, f32>>, k1: u32, k2: u32, v: f32) {
    data.entry(k1)
        .and_modify(|x| {
            x.entry(k2).and_modify(|x| *x = v).or_insert(v);
        })
        .or_insert_with(|| HashMap::from_iter([(k2, v)]));
}

pub fn phase_communities(
    data: HashMap<u32, HashMap<u32, f32>>,
    ref_weight: Option<HashMap<u32, f32>>,
) -> Vec<u32> {
    fn stat_ref_weight(ref_weight: &HashMap<u32, f32>, nodes: &HashSet<u32>) -> (i32, NotNan<f32>) {
        let mut count = 0;
        let mut weight = 0.;
        for node in nodes {
            if let Some(v) = ref_weight.get(node) {
                if *v > 0. {
                    count += 1;
                }else if *v < 0. {
                    count -= 1;
                }
                weight += v;
            }
        }
        (count, NotNan::new(weight).unwrap()) //safely
    }

    let lv = Louvain::new(data);
    let (data, mut communities) = lv.execute();

    if let Some(ref ref_weight) = ref_weight {
        // sort communities by the count and sum of weights between each community and reference from large to small,
        // this solution may need to be optimized
        communities.sort_by_cached_key(|x| Reverse(stat_ref_weight(&ref_weight, &x.nodes)));
    } else {
        // sort communities by weight from large to small
        communities.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap());
    }

    // let mut valid_ids = HashSet::default();
    let mut invalid_ids = HashSet::default();
    for (p, community) in communities.iter().enumerate() {
        if invalid_ids.contains(&community.id) {
            continue;
        }
        // valid_ids.insert(community.id);
        if let Some(id_vs) = data.get(&community.id) {
            for community_check in &communities[p + 1..] {
                if invalid_ids.contains(&community_check.id) {
                    continue;
                }
                if id_vs.contains_key(&community_check.id) {
                    invalid_ids.insert(community_check.id);
                }
            }
        }
    }

    // for c in &communities {
    //     println!("{} {} {} {} {:?} {:?}", c.id, !invalid_ids.contains(&c.id), c.weight, c.nodes.len(), 
    //         stat_ref_weight(&ref_weight.as_ref().unwrap(), &c.nodes), c.nodes);
    // }
    // panic!("");

    //merge nodes from the same phases
    let mut invalid_nodes = Vec::new();
    for community in communities
        .into_iter()
        .filter(|x| invalid_ids.contains(&x.id))
    {
        invalid_nodes.extend(community.nodes.into_iter());
    }
    invalid_nodes
}

fn init_graph(path: &str) -> HashMap<u32, HashMap<u32, f32>> {
    let mut data = new_data();
    let f = File::open(path).unwrap();
    for line in BufReader::new(f).lines() {
        let line = line.unwrap();
        if line == "\n" {
            continue;
        }
        let lines: Vec<&str> = line.split_ascii_whitespace().collect();
        insert_data(&mut data, u32::from_str(lines[1]).unwrap(), u32::from_str(lines[2]).unwrap(), f32::from_str(lines[3]).unwrap());
    }
    data
}

fn check_communities_are_correct(
    path: &str,
    data: &HashMap<u32, HashMap<u32, f32>>,
    communities: &[Node],
) {
    fn check_nodes_weight(data: &HashMap<u32, HashMap<u32, f32>>, node: &Node) {
        let mut w = 0.;
        for n1 in &node.nodes {
            for n2 in &node.nodes {
                w += data
                    .get(n1)
                    .map_or_else(|| 0., |v| *v.get(n2).unwrap_or(&0.));
            }
        }
        assert_eq!(w, node.weight * 2., "faield check nodes weight");
    }

    let raw_data = init_graph(path);
    for node1 in communities {
        check_nodes_weight(&raw_data, node1);
        for node2 in communities.iter().filter(|x| x.id != node1.id) {
            let mut w = 0.;
            for n1 in &node1.nodes {
                for n2 in &node2.nodes {
                    w += raw_data
                        .get(n1)
                        .map_or_else(|| 0., |v| *v.get(n2).unwrap_or(&0.));
                }
            }
            if w != 0. {
                assert_eq!(
                    w, data[&node1.id][&node2.id],
                    "faield check communities weight"
                );
            }
        }
    }
}

fn main() {
    let path: String = args().nth(1).unwrap();
    let lv = Louvain::new(init_graph(&path));
    let (data, communities) = lv.execute();
    check_communities_are_correct(&path, &data, &communities);
    for (p, c) in communities.iter().enumerate() {
        println!("{} {} {} {:?}", p, c.weight, c.nodes.len(), c.nodes);
    }
    let invalid_ids = phase_communities(init_graph(&path), Some(HashMap::default()));
    println!("{:?}", invalid_ids);
}
