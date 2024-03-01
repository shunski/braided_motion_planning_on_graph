use std::{collections::VecDeque, fmt::Debug, usize};

use augmented_graph::AugmentedGraph;
use topo_spaces::graph::RawSimpleGraph;
use util::SortedArray;

use std::ops;


pub mod graph_collection;
pub mod graphics;
pub mod search;
pub mod operators;
pub mod augmented_graph;
pub mod util;


#[derive(Clone, Copy, Debug)]
// This struct represents a cell in a discrete morse complex built upon the cube complex.
// 'N' is the number of points (or robots)
pub struct MorseCube<const N: usize> {
    // edges and vertices are stored in a separated variable. This implementation allows users to access the
    // edges (or vertices) quickly, even though the use of the memory is not the most efficient.

    // 'edges' is the set of edges in the cube.
    // The i-th element '[v_i, w_i]' of 'edges' represents the directed, critical motion of a point on vertex 'v_i' to 'w_i'.
    edges: util::SortedArray<N, [usize; 2]>,
    
    // 'vertices' is the set of vertices in the cube.
    vertices: util::SortedArray<N, usize>,

    // Note that, at any given time, the struct members must satisfy:
    //     (1) The union of 'edges' and 'vertices' must be a set of distinct points, and
    //     (2) The sum of 'self.edges.len()' and 'self.vertices.len()' is 'N'.
}


impl<const N: usize> MorseCube<N> {
    #[inline(always)]
    // 'fn new_unchecked' creates 'MorseCube' without performing any check on the input arguments.
    // This function reqiuires that;
    //     (1) elements of 'edges' and 'vertices' be distinct, and that
    //     (2) the sum of the size of 'edges' and 'vertices' be 'N'.
    // However, again, this method does not make sure these conditions are met.
    pub fn new_unchecked(edges: util::SortedArray<N, [usize; 2]>, vertices: util::SortedArray<N, usize>) -> Self {
        // elements of 'edges' and 'vertices' must be distinct, but we do not check this in this function.
        debug_assert!(
            {
                let mut v: Vec<_> = edges.iter().flatten().chain(vertices.iter()).collect();
                v.sort();
                v.iter().zip(v.iter().skip(1)).all(|(x,y)| x != y)
            }
        );
        // the size of 'edges' and 'vertices' must add up to 'N', but we do not check them in this function.
        debug_assert!( edges.len() + vertices.len() == N );
        
        Self { edges, vertices }
    }

    #[inline]
    // 'fn new_checked' creates 'MorseCube' from edges and vertices.
    // This function reqiuires that;
    //     (1) elements of 'edges' and 'vertices' be distinct, and that
    //     (2) the sum of the size of 'edges' and 'vertices' be 'N'.
    // These conditions will be checked internally. If the caller is 100% sure that the conditions are met and would like 
    // to improve the performance by skipping these checks, 'fn new_unchecked' can be used instead.
    pub fn new_checked(edges: util::SortedArray<N, [usize; 2]>, vertices: util::SortedArray<N, usize>) -> Self {
        // check that the elements of 'edges' and 'vertices' must be distinct. ToDo: This can be more efficient.
        assert!(
            {
                let mut v: Vec<_> = edges.iter().flatten().chain(vertices.iter()).collect();
                v.sort();
                v.iter().zip(v.iter().skip(1)).all(|(x,y)| x != y)
            }
        );
        // check that the size of 'edges' and 'vertices' must add up to 'N'.
        assert!( edges.len() + vertices.len() == N );
        
        Self { edges, vertices }
    }


    // 'fn flows_to_critical_one_cell' returns true iff 'self' flows to a critical one cell in the graph.
    pub fn flows_to_critical_one_cell(&self, graph: &RawSimpleGraph) -> bool {
        match self.edges.len() {
            0 => {
                // Any cube of dimension 0 does not flow to any critical 1-cell.
                false
            },

            1 => {
                let [initial, terminal] = self.edges[0];
                let idx1 = self.vertices.binary_search(&initial).expect_err("'edges' and 'vertices' are not distinct.");
                let idx2 = self.vertices.binary_search(&terminal).expect_err("'edges' and 'vertices' are not distinct.");
                
                // 'idx1==idx2' implies that the edge is order-respecting.
                idx1!=idx2
            },

            2 => {
                panic!("The case where the dimension is 2 is not supported yet.")
            },

            3.. => {
                // Any cube of dimension greater than 2 does not flow to any critical 1-cell.
                false
            },
        }
    }
    
    pub fn get_edge_path(&self, graph: &RawSimpleGraph) -> CubicPath<N> {
        let (critical_motion, [start, end]) = self.get_critical_path_and_its_ends();

        let mut path = VecDeque::from([critical_motion]);
        
        get_edge_path_recc::<N, false>(&start, &mut path, graph);
        get_edge_path_recc::<N, true>(&end, &mut path, graph);

        let path = CubicPath{
            path: path.into(),
            start,
            end
        };

        path
    }


    // 'fn get_critical_path_and_its_ends' returns the motion of the critical 1-cell as an 'ElementaryCubicPath' and its ends as "[SortedArray<N, usize>; 2]".
    // This method panics if 'self' is not a 1-dimensional cube. (Note that the dimension of 'self' is 'self.edges.len()')
    fn get_critical_path_and_its_ends(&self) -> (ElementaryCubicPath<N>, [SortedArray<N, usize>; 2]) {
        assert!(
            self.edges.len()==1, 
            "a cubic path can be created only if the cube is 1-dimensional, but it is {} dimensional. cube={self:?}", 
            self.edges.len()
        );

        let edge = self.edges[0];

        let mut start = self.vertices;
        start.insert(edge[0]);

        let mut end = self.vertices;
        end.insert(edge[1]);

        let out = {
            let mut out = ElementaryCubicPath::identity();
            out.data.push(edge);
            out
        };

        (out, [start, end])
    }
}


// 'ElementaryCubicPath' represents the one-step of the cubic path.
// This struct is PRIVATE to the module. Hence any method of this struct will not panic 
// even if the input argument is not meeting its condition. 
#[derive(Clone, Copy)]
struct ElementaryCubicPath<const N: usize> {
    // 'data' is a cube that represents a single (directed) step in the configuration space.
    // All the 2*'N' elements must be distinct.
    data: util::SortedArray<N,[usize;2]>
}

impl<const N: usize> Debug for ElementaryCubicPath<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data.fmt(f)
    }
}


impl<const N: usize> ElementaryCubicPath<N> {
    // 'fn identity' creates the identity path, that is, the path that does not move the element.
    fn identity() -> Self {
        Self{data: util::SortedArray::new() }
    }

    // 'fn is_identity' returns 'true' iff the path is the same as 'identity()'.
    fn is_identity(&self) -> bool {
        self.data.len()==0
    }


    // 'fn new_unchecked' creates the new 'ElementaryCubicPath' object from the input argument without performing any check on it.
    // The input argument 'v' must satisfy the requirements to be the member 'data' (See the struct definition).
    // However, this method does not check them.
    #[inline(always)]
    fn new_unchecked(v: util::SortedArray<N, [usize; 2]>) -> Self {
        // The input argument has to be sorted by the first element, by the order of usize, but we do not check this in this function.
        debug_assert!(
            v.iter().zip(v.iter().skip(1)).all(|([a,_],[b,_])| a < b),
            "The input argument has to be sorted by the first element, by the order of usize, but it is not;  {:?}",
            v
        );

        // The input argument has to have distinct elements, but we do not check this in this function.
        debug_assert!(
            {
                let mut elements: Vec<_> = v.iter().flatten().collect();
                elements.sort();
                elements.iter().zip(elements.iter().skip(1)).all(|(a,b)| a != b )
            },
            "The input argument has to have distinct elements, but it is not;  {:?}",
            v
        );


        ElementaryCubicPath{ data: v }
    }


    // 'fn new_unchecked' creates the new 'ElementaryCubicPath' object from the input argument without performing any check on it.
    // The input argument 'v' must satisfy the requirements to be the member 'data' (See the struct definition).
    // This method checks whether these conditions are met and panics if not. If the caller is 100% sure that these conditions are 
    // met and would like to improve the performance, then 'new_unchecked' might be used instead.
    fn new_checked(v: util::SortedArray<N, [usize; 2]>) -> Self {
        // The input argument has to be sorted by the first element, by the order of usize.
        assert!(
            v.iter().zip(v.iter().skip(1)).all(|([a,_],[b,_])| a < b),
            "The input argument has to be sorted by the first element, by the order of usize, but it is not;  {:?}",
            v
        );

        // The input argument has to have distinct elements.
        assert!(
            {
                let mut elements: Vec<_> = v.iter().flatten().collect();
                elements.sort();
                elements.iter().zip(elements.iter().skip(1)).all(|(a,b)| a != b )
            },
            "The input argument has to have distinct elements, but it is not;  {:?}",
            v
        );


        ElementaryCubicPath{ data: v }
    }


    // 'fn act_unchecked(&mut p)' takes in the mutable reference of distinct (but not necesarily sorted) points and moves 'p' 
    // by the action of the path 'self'. This method requires that;
    //     (1) 'p' be a set of DISTINCT points, and that
    //     (2) 'self[..][0]' be contained in 'p'.
    // However, it does not internally check these conditions. 
    fn act_unchecked(&self, p: &mut [usize; N]) {
        // 'self[..][0]' must be sorted by the order of 'usize'.
        debug_assert!(
            self.data.iter().map(|[a,_]| a)
                .zip(self.data.iter().map(|[a,_]| a))
                .all(|(x,y)| x < y),
            "'self.data[..][0]' must be sorted by the order of 'usize', but it is {:?}",
            self.data
        );

        // 'p' must be an array of (unsorted) distinct points
        debug_assert!({
            let mut q = p.clone();
            q.sort();
            q.iter().zip( q.iter().skip(1) ).all(|(&a, &b)| a < b )
        }, "'p' must be an array of (unsorted) distinct points, but p={p:?}");

        // in order for 'self' to be able to act on 'p', start of 'self' must be contained in 'p'
        debug_assert!(
            self.data.iter().all(|[s,_]| p.iter().any(|x| s==x ) ),
            "self={self:?}, p={p:?}"
        );

        // send each vertex to the goal
        for x in p {
            if let Ok(idx) = self.data.binary_search_by_key(x, |&[a, _]| a ) {
                *x = self.data[idx][1]
            };
        }
    }

    // 'fn act_setwise_unchecked()' takes in the mutable reference 'p' of distinct sorted points and moves it 
    // by the action of the path 'self' while keeping it sorted. This method requires that;
    //     (1) 'p' be a set of DISTINCT points, and that
    //     (2) 'self[..][0]' be contained in 'p'.
    // However, it does not internally check these conditions. 
    fn act_setwise_unchecked(&self, p: &mut SortedArray<N, usize>) {
        // 'p' must be an array of distinct points, but we do not check in this function.
        debug_assert!(
            p.iter().zip( p.iter().skip(1) ).all(|(&a, &b)| a < b ),
            "'p' must be an array of distinct points, but it is not: p={p:?}."
        );

        // In order for 'self' to be able to act on 'p', start of 'self' must be contained in 'p', but we do not check in this function.
        debug_assert!(
            self.data.iter().all(|[s,_]| p.binary_search(&s).is_ok() ),
            "self={self:?}, p={p:?}"
        );


        for [s,g] in self.data.iter() {
            let s_idx = p.binary_search(&s).ok().unwrap();
            p.modify(s_idx, *g);
        }
    }
}



// 'CubicPath' represents an 
#[derive(Clone)]
pub struct CubicPath<const N: usize> {
    // 'path' represents the a path from 'start' to 'end'.
    // This must be a "path" from 'start' to 'end', that is, 'self.act([usize; N]::from(self.start))' is setwise the same as 'end'.
    path: Vec<ElementaryCubicPath<N>>,
    
    // 'start' is the start configuration of the path.
    // It must have 'N' distinct elements.
    start: SortedArray<N, usize>,
    
    // 'end' is the end configuration of the path. 
    // It must have 'N' distinct elements.
    end: SortedArray<N, usize>,
    
    // 'end' is set-wise the same as 'self.act_on(&mut [usize, N]::try_from( start_usize ))'.
}

impl<const N: usize> ops::Deref for CubicPath<N> {
    type Target = Vec<ElementaryCubicPath<N>>;
    fn deref(&self) -> &Self::Target {
        &self.path
    }
}

impl<const N: usize> Debug for CubicPath<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "start = {:?}\n", self.start)?;
        write!(f, "end = {:?}\n", self.end)?;
        for p in self.iter() {
            write!(f, "{p:?}\n")?;
        }
        write!(f, "")
    }
}

impl<const N: usize> CubicPath<N> {
    // 'fn new_unchecked' creates a new "CubicPath" object from the two input arguments without checking them.
    // The input arguments 'start' and 'end' must be the sorted array of 'N' distinct elements, and the other input argument
    // 'path' must be well-defined as a path from 'start' to 'end'.
    // but function does not check these conditions.
    #[inline(always)]
    fn new_unchecked(path: Vec<ElementaryCubicPath<N>>, start: SortedArray<N, usize>, end: SortedArray<N, usize>) -> Self {
        // 'start' must have 'N' elements, but this method does not check it.
        debug_assert!(
            start.len() == N,
            "'start' must have 'N' elements, but it has only {} elements. 'start' = {:?}",
            start.len(),
            start
        );

        // The elements of 'start' must be distinct, but this method does not check it.
        debug_assert!(
            start.iter().zip(start.iter().skip(1)).all(|(a,b)| a!=b ),
            "The elements of 'start' must be distinct, but they are not: 'start' = {:?}",
            start
        );

        // 'end' must have 'N' elements, but this method does not check it.
        debug_assert!(
            end.len() == N,
            "'end' must have 'N' elements, but it has only {} elements. 'end' = {:?}",
            end.len(),
            end
        );

        // The elements of 'end' must be distinct, but this method does not check it.
        debug_assert!(
            end.iter().zip(end.iter().skip(1)).all(|(a,b)| a!=b ),
            "The elements of 'end' must be distinct, but they are not: 'end' = {:?}",
            end
        );

        // 'path' must be well-defined as a path from 'start' to 'end', but this method does not check it.
        debug_assert!(
            {
                let mut start: [usize; N] = start.as_array();
                for p in path { p.act_unchecked(&mut start); }
                start.sort();
                SortedArray::<N, usize>::from(start) == end
            }
        );

        Self { path, start, end }
    }


    // 'fn new_checked' creates a new "CubicPath" object from the two input arguments.
    // The input arguments 'start' and 'end' must be the sorted array of 'N' distinct elements, 
    // and this function panics if these conditions are not met.
    #[inline(always)]
    fn new_checked(path: Vec<ElementaryCubicPath<N>>, start: SortedArray<N, usize>, end: SortedArray<N, usize>) -> Self {
        // check that 'start' has 'N' elements.
        assert!(
            start.len() == N,
            "'start' must have 'N' elements, but it has only {} elements. 'start' = {:?}",
            start.len(),
            start
        );

        // check that the elements of 'start' are distinct.
        assert!(
            start.iter().zip(start.iter().skip(1)).all(|(a,b)| a!=b ),
            "The elements of 'start' must be distinct, but they are not: 'start' = {:?}",
            start
        );

        // check that 'end' has 'N' elements.
        assert!(
            end.len() == N,
            "'end' must have 'N' elements, but it has only {} elements. 'end' = {:?}",
            end.len(),
            end
        );

        // check that the elements of 'end' are distinct.
        assert!(
            end.iter().zip(end.iter().skip(1)).all(|(a,b)| a!=b ),
            "The elements of 'end' must be distinct, but they are not: 'end' = {:?}",
            end
        );


        // check that 'path' is well-defined as a path from 'start' to 'end'.
        debug_assert!(
            {
                let mut start: [usize; N] = start.as_array();
                for p in path { p.act_unchecked(&mut start); }
                start.sort();
                SortedArray::<N, usize>::from(start) == end
            }
        );

        Self { path, start, end }
    }


    // 'fn get_trivial_path_at' returns a trivial path at the input 'p', that is, path of length 0 that starts and ends at 'p'.
    // This method requires that 'p' be the set of 'N' distinct points and panics if this is not satisfied. 
    #[inline(always)]
    pub fn get_trivial_path_at(p: SortedArray<N, usize>) -> Self {
        // check that 'p' has 'N' elements.
        assert!(
            p.len() == N,
            "'p' must have 'N' elements, but it has only {} elements. 'p' = {:?}",
            p.len(),
            p
        );

        // check that the elements of 'p' are distinct.
        assert!(
            p.iter().zip(p.iter().skip(1)).all(|(a,b)| a!=b ),
            "The elements of 'p' must be distinct, but they are not: 'p' = {:?}",
            p
        );
        
        // 'path' will be the empty vector. This trivially satisfies the requirement for 'path'.
        Self { path: Vec::new(), start: p, end: p }
    }


    // 'fn act' takes in a configuration of 'N' points and computes the action of the path on the configuration.
    // This function requires that the input argument be a set of distinct (but not necessarily sorted) points.
    pub fn act(&self, p: &mut [usize; N]) {
        // check that the input argument be a set of distinct (but not necessarily sorted) points.
        assert!( {
                let mut p = p;
                p.sort();
                p.iter().zip(p.iter().skip(1)).all(|(a,b)| a < b)
            },
            "The input 'p' must be a set of distinct points, but they are not: p={:?}",
            p 
        );

        for e in self.iter() {
            // Note that having checked the condition above, the requirement for "fn ElementaryCubicPath::act_unchecked" 
            // is satisfied at each iteration.
            e.act_unchecked(p);
        }
    }

    // 'fn act' takes in a configuration of 'N' points and computes the action of the path on the configuration without performing checks.
    // This function requires that the input argument be a set of distinct (but not necessarily sorted) points, but it does not check the condition.
    pub fn act_unchecked(&self, p: &mut [usize; N]) {
        // the input argument must be a set of distinct (but not necessarily sorted) points, we do not check this.
        debug_assert!( {
                let mut p = p;
                p.sort();
                p.iter().zip(p.iter().skip(1)).all(|(a,b)| a < b)
            },
            "The input 'p' must be a set of distinct points, but they are not: p={:?}",
            p 
        );

        for e in self.iter() {
            // Note that having checked the condition above, the requirement for "fn ElementaryCubicPath::act_unchecked" 
            // is satisfied at each iteration.
            e.act_unchecked(p);
        }
    }

    pub fn composed_with(mut self, mut other: Self) -> Self {

        // 'self.end' must coinside with 'other.start'.
        debug_assert!(
            self.end == other.start,
            "cannot compose the two paths. path1 has end = {:?}, \nbut path2 has start = {:?}.",
            self.end, other.start
        );

        self.path.append(&mut other.path);


        self.reduce_to_geodesic()
    }


    // 'fn let_end_flow_to' computes the path from 'self.end' to 'new_end' and compose it to 'self', given the certain conditions.
    // This method guerantees that the resulting path be a geodesic if 'self' originally is.
    // This method requires that: 
    //     (1) 'new_end' have 'N' elements, that
    //     (2) the elements of 'new_end' be distinct, and that
    //     (3) each 'self.end[i]' is a child of 'new_end[i]'.
    // Internal checks on these will be performed in the method. If the user is 100% sure that these conditions are met,
    // then 'let_end_flow_to_unchecked' can be used instead.
    pub fn let_end_flow_to(&mut self, new_end: &SortedArray<N, usize>, graph: &AugmentedGraph) {
        // check that 'new_end' has 'N' elements.
        assert!(
            new_end.len() == N,
            "The input argument 'new_end' must have 'N' distinct elements, but is does not. new_end={new_end:?}"
        );

        
        // check that the elements of 'new_end' be distinct.
        assert!(
            new_end.iter().zip(new_end.iter().skip(1)).all(|(a,b)| a<b ),
            "input argument 'new_end' must be an array of distinct elements, but it is not: new_end={new_end:?}."
        );

        // check that each 'self.end[i]' is a child of 'new_end[i]'.
        assert!(
            self.end.iter().zip(new_end.iter()).all(|(&a,&b)| graph.is_first_child_of_second(a, b) )
        );


        // Now call 'fn let_end_flow_to_unchecked', which is the unchecked version of this funtion.
        // Because we have explicitly checked the condtions, it is safe to call the unchecked version.
        unsafe{
            self.let_end_flow_to_unchecked(new_end, graph)
        }
    }


    // 'fn let_end_flow_to' computes the path from 'self.end' to 'new_end' and compose it to 'self', given the certain conditions.
    // This method guerantees that the resulting path be a geodesic if 'self' originally is.
    // This method requires that: 
    //     (1) 'new_end' have 'N' elements, that
    //     (2) elements of 'new_end' be distinct, and that
    //     (3) each 'self.end[i]' is a child of 'new_end[i]'.
    // Internal checks on these conditions will not be performed in the method, and because some unsafe code depends on these conditions,
    // this function is unsafe.
    pub unsafe fn let_end_flow_to_unchecked(&mut self, new_end: &SortedArray<N, usize>, graph: &AugmentedGraph) {
        // The input 'new_end' must have 'N' elements, we do check in this funtion.
        debug_assert!(
            new_end.len() == N,
            "The input argument 'new_end' must have 'N' distinct elements, but is does not. new_end={new_end:?}"
        );
        
        // the elements of 'new_end' must be distinct, but we do not check in this function.
        debug_assert!(
            new_end.iter().zip(new_end.iter().skip(1)).all(|(a,b)| a<b ),
            "input argument 'new_end' must be an array of distinct elements, but it is not: new_end={new_end:?}."
        );

        // each 'self.end[i]' must be a child of 'new_end[i]', but we do not check in this function.
        debug_assert!(
            self.end.iter().zip(new_end.iter()).all(|(&a,&b)| graph.is_first_child_of_second(a, b) )
        );


        // The following while-loop succesively updates 'self.end' and 'self.end' will be 'new_end' at last.
        while &self.end != new_end {
            let mut elementary_path = SortedArray::<N, [usize; 2]>::new();

            // 'next_end' will be 'curr_end' in the next iteration. Its values will be filled in during the following for-loop (one value in each iteration),
            // and it will thus have length 'N' at the end of the present block. All the unsafe block that contains 'push_unchecked()' in the for-loop is
            // safe, because it will be called exactly 'N' times and because 'curr_end' is sorted.
            let mut next_end = SortedArray::new();


            // iterating over points.
            for (&s, &g) in self.end.iter().zip(new_end.iter()) {
                // 's' must flow to 'g' if they are not the same.
                if s == g {
                    unsafe{ next_end.push_unchecked(s) }; // see the intialization of 'next_end' about the safety of this block.
                    continue;
                }

                // find the vertex that 's' goes to next.
                // Note that 's-1' below does not panic, since 0 <= 'g' < 's'.
                let s_next = if graph.maximal_tree_contains([s-1, s]) {
                    // if 's-1' is connected by an edge in the maximal tree, then flow to it.
                    // Note that the edge ['s-1', 's'] must be order-respecting in this case. Hence, we do not need to do any check on this.
                    s-1
                } else {
                    // otherwise, 's' will flow to some essential vertex.
                    // Note that the following 'unwrap()' does not panic since 's' is adjacent to an essentail vertex smaller than 's'
                    // and since 's' is strictly positive.
                    let s_next = *graph.essential_vertices.iter().find(|&&x| graph.maximal_tree_contains([x, s]) ).unwrap();

                    // flow to 's_next' if the edge ['s', 's_next'] is order-respecting, or 's_next' occupied by another point of 'curr_end'. 
                    // Otherwise, skip it because we have to wait.
                    let result_i = self.end.binary_search(&s_next);
                    let j = self.end.binary_search(&s).err().unwrap(); // This should not panic, because 's' is in 'curr_end'.
                    if result_i.is_err() && result_i.err().unwrap() != j {
                        // This means that the edge ['s', 's_next'] is order-respecting, and 's_next' is not occupied because 'result_i' is not 'Ok(_)'. 
                        // Thus flow through the edge.
                        s_next
                    } else {
                        // Otherwise, update 'new_end' and skip.
                        unsafe{ next_end.push_unchecked(s) }; // see the intialization of 'next_end' about the safety of this block.
                        continue;
                    }
                };

                // update 'next_end'
                unsafe{ next_end.push_unchecked(s_next) }; // see the intialization of 'next_end' about the safety of this block.

                // add the step to the 'elementary_path'.
                // The following unsafe block is safe, because
                //     (1) this for-loop iterates exactly 'N' times (so this block is called at most 'N' times), and because
                //     (2) this '[s,s_next]' is always greater than the largest element as 'curr_end' is sorted.
                unsafe{ elementary_path.push_unchecked([s,s_next])}; 
            }

            self.path.push( ElementaryCubicPath::new_unchecked(elementary_path) );

            // Note that 'next_end' must have exactly 'N' elements at this point.
            debug_assert!(
                next_end.len() == N,
                "'next_end.len()' must be 'N'={N}, but it is {}",
                next_end.len()
            );

            // update 'curr_end'
            self.end = next_end;
        }
    }

    // 'fn let_end_reverse_flow_to' computes the path from 'self.end' to 'new_end' and composes it to 'self', given certain conditions are met.
    // There is the unique path determined by the discrete flow from 'new_end' to 'self.end', and the this method computes the inverse of it.
    // Note that computing the reverse flow is both simpler and more efficient than computing the flow because points do not need wait for other points
    // smaller than them (c.f. 'fn let_end_flow_to').
    // Furthermore, this method guerantees that the resulting path be a geodesic if 'self' originally is.
    // This method requires that: 
    //     (1) 'new_end' have 'N' elements, that
    //     (2) the elements of 'new_end' be distinct, and that
    //     (3) each 'new_end[i]' is a child of 'self.end[i]'.
    // Internal checks on these will be performed in the method. If the user is 100% sure that these conditions are met,
    // then 'let_end_reverse_flow_to_unchecked' can be used instead.
    pub fn let_end_reverse_flow_to(&mut self, new_end: &SortedArray<N, usize>, graph: &AugmentedGraph) {
        // check that 'new_end' has 'N' elements.
        assert!(
            new_end.len() == N,
            "The input argument 'new_end' must have 'N' distinct elements, but is does not. new_end={new_end:?}"
        );

        // check that the elements of 'new_end' are distinct.
        assert!(
            new_end.iter().zip(new_end.iter().skip(1)).all(|(a,b)| a<b ),
            "The input argument 'new_end' must be an array of distinct elements, but it is not: new_end={new_end:?}."
        );

        // check that each 'new_end[i]' is a child of 'self.end[i]'.
        assert!(
            new_end.iter().zip(self.end.iter()).all(|(&a,&b)| graph.is_first_child_of_second(a, b) )
        );


        // Now call 'fn let_end_reverse_flow_to_unchecked', which is the unchecked version of this funtion.
        // Because we have explicitly checked the condtions, it is safe to call the unchecked version.
        unsafe{
            self.let_end_reverse_flow_to_unchecked(&new_end, graph);
        }
    }


    // 'fn let_end_reverse_flow_to_unchecked' computes the path from 'self.end' to 'new_end' and composes it to 'self', given certain conditions are met.
    // There is the unique path determined by the discrete flow from 'new_end' to 'self.end', and the this method computes the inverse of it.
    // Note that computing the reverse flow is both simpler and more efficient than computing the flow because points do not need wait for other points
    // smaller than them (c.f. 'fn let_end_flow_to').
    // Furthermore, this method guerantees that the resulting path be a geodesic if 'self' originally is.
    // This method requires that: 
    //     (1) 'new_end' have 'N' elements, that
    //     (2) the elements of 'new_end' be distinct, and that
    //     (3) each 'new_end[i]' is a child of 'self.end[i]'.
    // Internal checks on these conditions will not be performed in the method, and because some unsafe code depends on these conditions,
    // this function is unsafe.
    pub unsafe fn let_end_reverse_flow_to_unchecked(&mut self, new_end: &SortedArray<N, usize>, graph: &AugmentedGraph) {
        // 'new_end' must have 'N' elements, but we do not check this in this function.
        debug_assert!(
            new_end.len() == N,
            "The input argument 'new_end' must have 'N' distinct elements, but is does not. new_end={new_end:?}"
        );

        // The elements of 'new_end' must be distinct, but we do not check this in this function.
        debug_assert!(
            new_end.iter().zip(new_end.iter().skip(1)).all(|(a,b)| a<b ),
            "The input argument 'new_end' must be an array of distinct elements, but it is not: new_end={new_end:?}."
        );

        // Each 'new_end[i]' must be a child of 'self.end[i]', but we do not check this in this function.
        debug_assert!(
            new_end.iter().zip(self.end.iter()).all(|(&a,&b)| graph.is_first_child_of_second(a, b) )
        );


        // The following while-loop succesively updates 'self.end', and 'self.end' will be 'new_end' at last.
        while &self.end != new_end {
            let mut elementary_path = SortedArray::<N, [usize; 2]>::new();

            // 'next_end' will be 'curr_end' in the next iteration. Its values will be filled in during the following for-loop (one value in each iteration),
            // and it will thus have length 'N' at the end of the present block. All the unsafe block that contains 'push_unchecked()' in the for-loop is
            // safe, because it will be called exactly 'N' times and because 'curr_end' is sorted.
            let mut next_end = SortedArray::new();


            // Iterating over points.
            for (&s, &g) in self.end.iter().zip(new_end.iter()) {
                // 's' must flow to 'g' if they are not the same.
                if s == g {
                    // update 'next_end'
                    unsafe{ next_end.push_unchecked(s) }; // see the intialization of 'next_end' about the safety of this block.
                    continue;
                }

                // Find the vertex that 's' goes to next.
                let s_next = if graph.essential_vertices.binary_search(&s).is_ok() {
                    // If 's' is an essential vertex then choose the branch that contains 'g'.
                    *graph.next_vertices[s].iter().take_while(|next| &s >= next ).last().unwrap()
                } else {
                    // If 's' is not an essential vertex, then it can simply go to the vertex 's+1'.
                    s+1
                };
                
                // update 'next_end'
                unsafe{ next_end.push_unchecked(s_next) }; // see the intialization of 'next_end' about the safety of this block.

                // Add the step to the 'elementary_path'.
                // The following unsafe block is safe, because
                //     (1) this for-loop iterates exactly 'N' times (so this block is called at most 'N' times), and because
                //     (2) this '[s,s_next]' is always greater than the largest element as 'curr_end' is sorted.
                unsafe{ elementary_path.push_unchecked([s,s_next])}; 
            }

            self.path.push( ElementaryCubicPath::new_unchecked(elementary_path) );

            // Note that 'next_end' must have exactly 'N' elements at this point.
            debug_assert!(
                next_end.len() == N,
                "'next_end.len()' must be 'N'={N}, but it is {}",
                next_end.len()
            );

            // update 'curr_end'
            self.end = next_end;
        }
    }

    fn reduce_to_geodesic(self) -> Self {
        let (start, end) = (self.start, self.end);
        let mut out = Vec::new();

        for f in self.path {
            if out.is_empty() {
                out.push(f);
                continue;
            }

            let mut new_motion = ElementaryCubicPath::identity();
            for [v, w] in f.data.into_iter() {
                if let Some(non_commuting_idx) = (0..out.len()).rev().find(|&i| out[i].data.iter().any(|&[x, y]| v==x || v==y || w==x || w==y ) ) {
                    let non_commuting_motion = &mut out[non_commuting_idx];
                    if let Ok(idx) = non_commuting_motion.data.binary_search(&[w, v] ) {
                        // if 'non_commuting_motion' contains the inverse motion of [v,w], remove it from 'non_commuting_motion'
                        non_commuting_motion.data.remove(idx);
                        
                        // remove 'non_commuting_motion' if the motion is now the identity motion
                        if non_commuting_motion.is_identity() {
                            out.remove(non_commuting_idx);
                        }
                    } else {
                        // if 'non_commuting_motion' does not contain the inverse
                        if non_commuting_idx == out.len()-1 {
                            new_motion.data.push([v,w]);
                        } else {
                            out[non_commuting_idx+1].data.insert([v,w]);
                        }
                    }
                } else {
                    // if '[v,w]' commutes with every motion in 'out', then insert the motion at the front.
                    out[0].data.insert([v,w]);
                }
            }

            // add 'new_motion' to 'out' if 'new_motion' is not the identity motion
            if !new_motion.is_identity() {
                out.push(new_motion);
            }

        }

        Self{
            path: out,
            start,
            end
        }
    }


    // 'fn get_geodesic_between' computes the geodesic between the two (ordered) configurations 'p' and 'q' in the given graph.
    // This method requires that 'p' and 'q' have exactly 'N' elements and panics if the condition is not met.
    fn get_geodesic_between([p, q]: [&SortedArray<N, usize>; 2], graph: &AugmentedGraph) -> Self {
        // check that both 'p' and 'q' have 'N' elements.
        assert!(
            p.len() == N,
            "The input argumet 'p' must have 'N' elements, but it does not: p={p:?}"
        );
        assert!(
            q.len() == N,
            "The input argumet 'q' must have 'N' elements, but it does not: q={q:?}"
        );

        // 'meets_without_multiplicity' is the meets of points in 'p' and 'q' with mutiplicity reduced.
        let meets_without_multiplicity: SortedArray<N,_> = {
            // Collect the "meets". Note that 'meets' is already sorted because 'p' and 'q' are sorted.
            let mut meets: Vec<_> = p.iter().zip(q.iter()).map(|(&s, &g)| graph.meet_of(s, g) ).collect();
            

            // Reduce multiplicity of 'meets'.
            let mut multiplicity_count = 0;
            let mut prev_meet = meets.first().unwrap() + 1; // initial value of 'prev_meet' only needs to be different from 'meets.first().unwrap()'.
            for meet in &mut meets {
                if &prev_meet == meet {
                    // If 'meet' is the same as the previous 'meet' ('prev_meet'), then substract 'meet' by 'multiplicity_count'.
                    *meet -= multiplicity_count;

                    // update 'multiplicity_count'.
                    multiplicity_count += 1;
                } else {
                    // If 'meet' is not the same as the previous 'meet' ('prev_meet'), then no change will be made to 'meet'.

                    // reset 'multiplicity_count'.
                    multiplicity_count = 1;
                    
                    // update 'prev_meet'
                    prev_meet = *meet;
                }
            }

            // At this point, 'meets' is sorted and contains 'N' distinct elements.
            debug_assert!(meets.iter().zip(meets.iter().skip(1)).all(|(i,j)| i<j ));
            debug_assert!(meets.len()==N);

            // Construct the SortedArray. This does not panic because 'meets' has 'N' elements and is sorted.
            // The following "try_from()" does not panic because 'meets' is sorted.
            // The following "err().unwrap()" does not panic because 'meets' has exactly 'N' elements.
            SortedArray::try_from( meets ).err().unwrap()
        };

        // Compute the first half of the output path, that is, the path from 'p' to 'meets_without_multiplicity'.
        // The following unsafe block is safe, because 
        //     (1) 'meets_without_multiplicity' has 'N' elements, 
        //     (2) the elements of 'meets_without_multiplicity' are distinct (this is because it is free from multiplicity), and
        //     (3) each 'p[i]' is a child of 'meets_without_multiplicity[i]'. 
        let mut first_half = CubicPath::get_trivial_path_at(*p);
        unsafe {
            first_half.let_end_flow_to_unchecked(&meets_without_multiplicity, graph)
        };


        // Compute the first half of the output path, that is, the path from 'p' to 'meets_without_multiplicity'.
        // The following unsafe block is safe, because 
        //     (1) 'meets_without_multiplicity' has 'N' elements, 
        //     (2) the elements of 'meets_without_multiplicity' are distinct (this is because it is free from multiplicity), and
        //     (3) each 'q[i]' is a child of 'meets_without_multiplicity[i]'.
        let mut second_half = CubicPath::get_trivial_path_at(meets_without_multiplicity);
        unsafe {
            second_half.let_end_reverse_flow_to_unchecked(&meets_without_multiplicity, graph)
        };


        // Compose the two paths and return it.
        first_half.composed_with(second_half)

    }



    // 'fn extend_by_checked' extends the path 'self' by concatenating a unit path 'path' at the end.
    // This method requires that the start vertices of 'path' be contained in the end of 'self' and panics if the condition is not met.
    fn extend_by_checked(&mut self, path: ElementaryCubicPath<N>) {
        // check that the start vertices of 'path' be contained in the end of 'self'.
        assert!(
            path.data.iter().all(|[s,_]| self.end.binary_search(s).is_ok() ),
            "The start vertices of 'path' is not contained in the end of 'self'."
        );

        // Now push 'path' to 'self.path', and modify 'self.end'.
        self.path.push( path );
        path.act_setwise_unchecked( &mut self.end );
    }


    // 'fn extend_by_unchecked' extends the path 'self' by concatenating a unit path 'path' at the end.
    // This method requires that the start vertices of 'path' be contained in the end of 'self', but this function does not check the condition.
    fn extend_by_unchecked(&mut self, path: ElementaryCubicPath<N>) {
        // The start vertices of 'path' must be contained in the end of 'self'.
        debug_assert!(
            path.data.iter().all(|[s,_]| self.end.binary_search(s).is_ok() ),
            "The start vertices of 'path' is not contained in the end of 'self'."
        );

        // Now push 'path' to 'self.path', and modify 'self.end'.
        self.path.push( path );
        self.act_setwise_unchecked( &mut self.end );
    }
}

pub struct MorsePath<'a, const N: usize> {
    // 'start' is the start configuration of the path.
    // It always has 'N' distinct elements.
    start: SortedArray<N , usize>,

    // 'end' is the end configuration of the path.
    // It always has 'N' distinct elements.
    end: SortedArray<N , usize>,

    path: Vec<MorseCube<N>>,
    graph: &'a AugmentedGraph<'a>,
}

impl<'a, const N: usize> MorsePath<'a, N> {

    // 'fn new' creates a new "MorsePath" object from the two input arguments.
    // The input arguments must be a sorted array of 'N' distinct elements, and this function panics if these conditions are not met. 
    pub fn new(start: SortedArray<N, usize>, end: SortedArray<N, usize>, graph: &'a AugmentedGraph) -> Self {
        // check that 'start' has 'N' elements.
        assert!(
            start.len() == N,
            "'start' must have 'N' elements, but it has only {} elements. 'start' = {:?}",
            start.len(),
            start
        );

        // check that the elements of 'start' are distinct.
        assert!(
            start.iter().zip(start.iter().skip(1)).all(|(a,b)| a!=b ),
            "The elements of 'start' must be distinct, but they are not: 'start' = {:?}",
            start
        );

        // check that 'end' has 'N' elements.
        assert!(
            end.len() == N,
            "'end' must have 'N' elements, but it has only {} elements. 'end' = {:?}",
            end.len(),
            end
        );

        // check that the elements of 'end' are distinct.
        assert!(
            end.iter().zip(end.iter().skip(1)).all(|(a,b)| a!=b ),
            "The elements of 'end' must be distinct, but they are not: 'end' = {:?}",
            end
        );

        Self {
            start,
            end,
            path: Vec::new(),
            graph,
        }
    }

    pub fn add(&mut self, c: MorseCube<N>) {
        self.path.push(c);
    }


    // 'fn get_geodesic' computes the geodesic of the path described by the morse path 'self'.
    pub fn get_geodesic(&'a self) -> CubicPath<N> {
        self.path.iter()
            .fold(
                CubicPath::get_trivial_path_at(self.start) , 
                |accum, cube| {
                    let start = accum.end;
                    let (critical_path, [end, _]) = cube.get_critical_path_and_its_ends();
                    accum.composed_with(CubicPath::get_geodesic_between([&start, &end], self.graph));

                    accum.extend_by(critical_path);

                    // return
                    accum
                } 
            )
            .composed_with(other);
    }
}


#[cfg(test)]
mod morse_path_test {
    // 'N' is the number of robots
    const N: usize = 5;
    use crate::{ graph_collection, graphics, MorsePath, MorseCube, util::SortedArray };

    // #[test]
    fn morse_path() -> Result<(), Box<dyn std::error::Error>> {
        let (graph, embedding, name) = graph_collection::RawGraphCollection::Grid.get();

        let cube1 = {
            let edges = SortedArray::from([[70,105]]);
            let vertices = SortedArray::from([71,72,74,75]);
            MorseCube::<N>::new_unchecked(edges, vertices)
        };
        cube1.is_valid(&graph);

        let cube2 = {
            let edges = SortedArray::from([[88,123]]);
            let vertices = SortedArray::<N, usize>::from([89,90,92,93]);
            MorseCube::<N>::new_unchecked(edges, vertices)
        };
        cube2.is_valid(&graph);
        
        let start = [11, 26, 56, 67, 132];
        let end = [3, 17, 100, 110, 120];
        let morse_path = MorsePath{ start, end, path: vec![cube1, cube2], graph: &graph };

        let path = morse_path.get_geodesic();

        // draw the generated paths
        println!("drawing paths...");
        graphics::draw_edge_path::<N>(&path, &"morse_path_test", &name, &embedding, &graph)?;

        Ok(())
    }
}


// 'fn get_edge_path_recc' recursively computes the path that connects 'points' and the critical 0-cell.
// The resulting path is stored in the 'path'.
// The direction of the path should be specified by the caller by the 'FORWARD' value; if 'FORWARD' is true (resp. false), then the path 
// will be created in such a way that the critical 0-cell comes at first (resp. at last). 
fn get_edge_path_recc<const N: usize, const FORWARD: bool>(points: &SortedArray<N, usize>, path: &mut VecDeque<ElementaryCubicPath<N>>, graph: &RawSimpleGraph) {
    // The base case: If the path is created all the way to the critical 0-cell, then return.
    if points.iter().enumerate().all(|(i, &x)| i==x ) {
        return;
    }

    let (elementary_path, next_points) = {
        let mut out = util::SortedArray::new();
        let mut next_points = SortedArray::new();
        for &p in points.iter() {
            let next = graph.adjacent_vertices_iter(p)
                .take_while(|&v| v < p)
                .filter( |&v| graph.maximal_tree_contains([v, p]) )
                .next().unwrap_or(p);
            
            // check whether 'next' is occupied by some point from 'points' or not.
            // if it does, then 'p' cannot proceed to 'next'.
            let next = if points.binary_search(&next).is_ok() {
                p
            } else {
                // Now 'p' can proceed to 'next' if it is order-respecting.
                if next+1 == p {
                    next
                } else if points.iter().all(|q| !(next..p).contains(q) ) { // TODO: this can be more efficient
                    next
                } else {
                    // Otherwise 'p' cannnot proceed to 'next'.
                    p
                }
            };


            // If 'p' is moving to another vertex, record the motion.
            if p != next {
                if FORWARD {
                    out.insert([p, next]);
                } else {
                    out.insert([next, p]);
                };
            }

            // record 'next'
            next_points.push( next );
        }
        (ElementaryCubicPath{ data: out }, next_points)
    };

    // record the motion in 'path'.
    if FORWARD {
        path.push_back( elementary_path );
    } else {
        path.push_front( elementary_path );
    }


    // run the function recursively
    get_edge_path_recc::<N, FORWARD>(&next_points, path, graph)
}