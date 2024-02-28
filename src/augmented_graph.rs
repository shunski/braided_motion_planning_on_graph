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


    // 'fn meet_of' computes the "meet" of two vertices 'v1' and 'v2' in the chosen maximal tree.
    // The "meet" of two vertices is defined by as the smallest vertex contained in the gwodesic between the two vertices.
    pub fn meet_of(&self, v: usize, w: usize) -> usize {

        // If one vertex is the child of the other, then the meet is smaller vertex.
        if v == w {
            return v;
        } else if v < w && self.is_first_child_of_second(w, v) {
            // If 'w' is the child of 'v'
            return v;
        } else if w < v && self.is_first_child_of_second(v, w) {
            return w;
        }

        // Now we can assume that the meet is an essential vertex.

        // convert 'v' and 'w' to the supremum of essential vertices less than them.
        let mut v = match self.essential_vertices.binary_search(&v) {
            Ok(idx) => self.essential_vertices[idx],
            Err(idx) => self.essential_vertices[idx-1],
            // This does not panic because '0' is an essential vertex
        };

        let mut w = match self.essential_vertices.binary_search(&w) {
            Ok(idx) => self.essential_vertices[idx],
            Err(idx) => self.essential_vertices[idx-1],
            // This does not panic because '0' is an essential vertex
        };

        // Now find the meet by tracing down the parents.
        while v != w {
            let greater_point = if v > w {
                &mut v
            } else {
                &mut w
            };

            *greater_point = self.parent[*greater_point];
        }

        v
    }


    // 'fn is_first_child_of_second' determines whether the first input ('v') is the child of the second input ('w').
    // A vertex 'v' is a "child" of a vertex 'w' if 'v'<'w' and 'v' is the smallest vertex contained in the geodesic between them in the tree.
    pub fn is_first_child_of_second(&self, mut v: usize, w: usize) -> bool {
        // if 'v' is a child of 'w', then 'v' > 'w'.
        if v <= w {
            return false;
        }

        // 'x' is the smallest essential vertex greater than or equal to 'w'.
        let x = match self.essential_vertices.binary_search(&w) {
            Ok(idx) => idx,
            Err(idx) => idx,
        };

        // If the geodesic between 'v' and 'w' does not contain any essential vertex, then 'v' is a child of 'w'. 
        if v < x {
            return true;
        }

        // Now 'v' >= max{ 'w', 'x' }.
        let parent = self.parent[x];

        // 'sup_of_children' is such that a vertex 'y' is contained in the range ('v'..'sup_of_children') iff 'y' is the child of 'v'.
        let sup_of_children = match self.next_vertices[parent].binary_search(&x) {

            // 'Ok(_)' means 'x' and its parent is adjacent (which rarelly happens because both are essential).
            // In this case, we have to choose 'sup_of_childs' to be 'self.next_vertices[parent][idx+1]'. If such choice is not possible,
            // then choose 'usize::MAX'
            Ok(idx) => if idx+1 < self.n_vertices() {
                self.next_vertices[parent][idx+1]
            } else {
                usize::MAX
            },

            // 'Err(idx)' means that 'self.next_vertices[parent][idx]' is the upper bound. If such choice is not possible,
            // then choose 'usize::MAX'.
            Err(idx) => if idx < self.n_vertices() {
                self.next_vertices[parent][idx]
            } else {
                usize::MAX
            },
        };

        (v..sup_of_children).contains(&w)
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
        // Just making sure that 'essential_vertices' is sorted...
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