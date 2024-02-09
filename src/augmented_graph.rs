use topo_spaces::graph::RawSimpleGraph;

fn get_next_essential_vertices(essential_vertex: usize, graph: &RawSimpleGraph) -> Vec<usize> {
    let mut out = Vec::new();
    for mut v in graph.adjacent_vertices_iter(essential_vertex)
        .skip(1)
        .filter(|&v| v > essential_vertex)
        .filter(|&v| graph.maximal_tree_contains( [essential_vertex, v] )) {
        while graph.degree_in_maximal_tree(v) == 2 {
            v += 1;
        }

        out.push(v);
    }
    out
}

fn get_next_vertices(essential_vertex: usize, graph: &RawSimpleGraph) -> Vec<usize> {
    graph.adjacent_vertices_iter(essential_vertex)
        .skip(1)
        .filter(|&v| v > essential_vertex)
        .filter(|&v| graph.maximal_tree_contains( [essential_vertex, v] ))
        .collect::<Vec<_>>()
}

struct UsizeIndexable<T>(Vec<(usize, T)>);

impl<T> std::ops::Index<usize> for UsizeIndexable<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[self.0.binary_search_by_key(&index, |&(i,_)| i ).unwrap()].1
    }
}

pub struct AugmentedGraph<'a> {
    // 'graph_ref' is the reference of the raw graph. The rest of the members of this struct is the information about
    // this raw graph.
    pub graph_ref: &'a RawSimpleGraph,


    //
    pub next_essential_vertices: UsizeIndexable<Vec<usize>>,
    

    // 
    pub next_vertices: UsizeIndexable<Vec<usize>>,


    // 'parent' stores all the "essential parent" of each essential vertex, that is, if 'v' is the index of 
    // an essential vertex in 'graph_ref', then 'parent[v]' is the maximum of those essential vertices that are
    // less than 'v' (Here, the order between the vertices are determined by the order of the indeces).
    // Note that this definition of parent is well-defined for all of the essential vertices except the vertex '0'.
    // Hence, 'parent[0]' will always panic.
    pub parent: UsizeIndexable<usize>,
    

    // 'essential_vertices' is the collection of all essential vertices,
    // whice are by definition, the vertices of degree not equal to 2 in maximal tree.
    pub essential_vertices: Vec<usize>, 
}

impl<'a> std::ops::Deref for AugmentedGraph<'a> {
    type Target = RawSimpleGraph;
    fn deref(&self) -> &Self::Target {
        &self.graph_ref
    }
}

impl<'a> AugmentedGraph<'a> {
    pub fn from(graph_ref: &'a RawSimpleGraph) -> Self {
        let essential_vertices = Self::get_essential_vertices(graph_ref);
        Self {
            graph_ref,
            next_essential_vertices: Self::get_next_essential_vertices_dictionary(graph_ref, &essential_vertices),
            next_vertices: Self::get_next_vertices_dictionary(graph_ref, &essential_vertices),
            parent: Self::get_parent(&essential_vertices),
            essential_vertices,
        }
    }

    fn get_next_essential_vertices_dictionary(graph: &RawSimpleGraph, essential_vertices: &[usize]) -> UsizeIndexable<Vec<usize>> {
        let out = essential_vertices.iter()
            .map(|&v| (v, get_next_essential_vertices(v, graph)))
            .collect::<Vec<_>>();

        UsizeIndexable(out)
    }
    
    
    fn get_next_vertices_dictionary(graph: &RawSimpleGraph, essential_vertices: &[usize]) -> UsizeIndexable<Vec<usize>> {
        let out = essential_vertices.iter()
            .map(|&v| (v, get_next_vertices(v, graph)))
            .collect::<Vec<_>>();

        UsizeIndexable(out)
    }
    
    fn get_essential_vertices(graph: &RawSimpleGraph) -> Vec<usize> {
        graph.vertex_iter()
            .map(|v| v.vertex() )
            .filter(|v| graph.degree_in_maximal_tree(*v) != 2 )
            .collect::<Vec<_>>()
    }

    fn get_parent(essential_vertices: &[usize]) -> UsizeIndexable<usize> {
        // Just making sure that the 'essential_vertices' are sorted...
        assert!(
            essential_vertices.iter().zip(&essential_vertices[1..]).all(|(v, w)| v < w ),
            "'essential_vertices' is not sorted."
        );

        let out = essential_vertices.iter()
            .skip(1) // skip the first vertex, which does not have a parent.
            .map(|&v| (v, *essential_vertices.iter().take_while(|&&w| w<v).last().unwrap()) )
            .collect::<Vec<_>>();
        UsizeIndexable(out)
    }
}