use std::{collections::VecDeque, fmt::Debug};

use topo_spaces::graph::RawSimpleGraph;

use std::ops::{Deref, DerefMut};
use std::{mem, usize};


pub mod graph_collection;
pub mod graphics;
pub mod search;
pub mod operators;


#[derive(Clone, Copy)]
pub struct SortedArray<const N: usize, T: Ord> {
    data: [T; N],
    len: usize,
}

impl<const N: usize, T: Ord + Copy> SortedArray<N, T> {}

impl<const N: usize, T: Ord + Copy> std::ops::Deref for SortedArray<N, T> {
    type Target = [T];
    fn deref(&self) -> &Self::Target {
        &self.data[..self.len]
    }
}

impl<const N: usize, T: Ord + Copy> std::ops::DerefMut for SortedArray<N, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data[..self.len]
    }
}

impl<const N: usize, T: Ord + Copy + Debug> Debug for SortedArray<N, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data[..self.len].fmt(f)
    }
}

impl<const N: usize, const M: usize, T: Ord + Copy + Debug> From<[T; M]> for SortedArray<N, T> {
    fn from(mut value: [T; M]) -> Self {
        assert!( M <= N, "cannot build an array of size {N} from an array of size {M}");
        
        value.sort();
        let mut out = Self::new();
        for x in value {
            out.add_without_sort(x);
        }
        out
    }
}


impl<const N: usize, T: Ord + Copy + Debug> TryFrom<Vec<T>> for SortedArray<N, T> {
    type Error = Self;
    // if 'value.len() > N', then this function creates an array of 'N' elements that contains first 'N' elements of 'value',
    // and enclose it in 'Err'
    fn try_from(mut value: Vec<T>) -> Result<Self, Self::Error> {
        value.sort();
        let mut out = Self::new();
        let is_error = value.len() > N;
        for x in value.into_iter().take(N) {
            out.add_without_sort(x);
        }
        if is_error {
            Ok(out)
        } else {
            Err(out)
        }
    }
}

impl<const N: usize, T: Ord + Copy + Debug> SortedArray<N,T> {
    pub fn new() -> Self {
        let data = unsafe {
            mem::MaybeUninit::uninit().assume_init()
        };
        Self { data, len: 0 }
    }

    #[inline]
    pub fn add_without_sort(&mut self, v: T) {
        assert!(self.len<N, "The array is full. array={self:?}, while trying to add {v:?}");

        unsafe{
            self.data.as_mut_ptr().add(self.len).write(v);
        }

        self.len+=1;
        debug_assert!(self.iter().zip(self.iter().skip(1)).all(|(x,y)| x<y ))
    }

    pub fn add(&mut self, v: T) {
        assert!(self.len<N, "The array is full. array={self:?}, while trying to add {v:?}");

        unsafe{
            self.data.as_mut_ptr().add(self.len).write(v);
        }
        self.len+=1;

        self.sort(); // ToDo: This can be more efficient
    }

    pub fn remove(&mut self, idx: usize) {
        assert!(idx < self.len, "index out of bounds: array is of length={}, but 'idx'={idx}", self.len);

        for i in idx..self.len-1 {
            self[i] = self[i+1]; 
        }

        self.len -= 1;
    }

    fn into_iter(self) -> impl Iterator<Item = T> { 
        self.data.into_iter().take(self.len)
    }
}


#[derive(Clone, Copy, Debug)]
// This struct represents a cell in a discrete morse complex built upon the cube complex.
// 'N' is the number of points (or robots)
pub struct MorseCube<const N: usize> {
    // edges and vertices are stored in a separated variable. This implementation allows users to access the
    // edges (or vertices) quickly, even though the memory space is not minimal.

    // 'edges' is the set of edges in the cube.
    edges: SortedArray<N, [usize; 2]>,
    
    // 'vertices' is the set of vertices in the cube.
    vertices: SortedArray<N, usize>,

    // It must be true that 'N = self.edges.len() + self.vertices.len()' at any given time.
}


impl<const N: usize> MorseCube<N> {
    #[inline]
    pub fn new_unchecked(edges: SortedArray<N, [usize; 2]>, vertices: SortedArray<N, usize>) -> Self {
        // the size of 'edges' and 'vertices' must add up to 'N', but we do not check them in this function.
        debug_assert!( edges.len() + vertices.len() == N );
        Self { edges, vertices }
    }

    // this function creates a morse cube only if the cell described by the input is critical. Otherwise it returns 'None'
    pub fn new_checked(edges: SortedArray<N, [usize; 2]>, vertices: SortedArray<N, usize>) -> Option<Self> {
        // the size of 'edges' and 'vertices' must add up to 'N', but we do not check them in this function.
        debug_assert!( edges.len() + vertices.len() == N );

        for mut edge in edges.iter().copied() {
            edge.sort();

            // The following closure is not necessary, but it can help speed up.
            if edge[1] - edge[0] == 1 {
                return None;
            };

            // is the edge is not disrespecting the order then return 'None'
            if vertices.iter().copied().all(|v| v < edge[0] || edge[1] < v ) {
                return None;
            }
        }

        // coming here means that the edge is critical.
        Some( Self { edges, vertices } )
    }

    pub fn is_valid(&self, graph: &RawSimpleGraph) -> bool {
        // check that the number of points is correct
        if self.edges.len() + self.vertices.len() != N {
            return false;
        }

        // check that the cell is contained in the graph
        self.vertices.iter().all( |&v| v < graph.n_vertices() )
        && self.edges.iter().all(|&[v, w]| graph.contains(v, w) )
    }
    
    pub fn get_edge_path(&self, graph: &RawSimpleGraph) -> CubicPath<N> {
        let (critical_motion, [start, end]) = ElementaryCubicPath::from_morse_cube(self);

        // 'start' and 'end' are already sorted.
        debug_assert!( 
            end.into_iter().zip(end.into_iter().skip(1)).all(|(a,b)| a < b ),
            "end = {end:?} \n is not ordered. implementation of 'ElementaryCubicPath::from_morse_cube' might be wrong."
        );
        debug_assert!( 
            start.into_iter().zip(start.into_iter().skip(1)).all(|(a,b)| a < b ),
            "start = {start:?} \n is not ordered. implementation of 'ElementaryCubicPath::from_morse_cube' might be wrong."
        );

        let mut path = VecDeque::from([critical_motion]);
        // panic!("start={start:?}");
        // panic!("end={end:?}");
        let start = get_edge_path_recc::<N, false>(&start, &mut path, graph);
        let end = get_edge_path_recc::<N, true>(&end, &mut path, graph);

        let path = CubicPath{
            path: path.into(),
            start,
            end
        };

        path
    }
}


// This struct represents the one-step of the cubic path.
// Any instance of this struct must satisfy the condition that the array 'self[..][0]' is sorted, that is, the start of elementary path is sorted. 
#[derive(Clone, Copy)]
pub struct ElementaryCubicPath<const N: usize>(SortedArray<N,[usize;2]>);

impl<const N: usize> Debug for ElementaryCubicPath<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl<const N: usize> Deref for ElementaryCubicPath<N>{
    type Target = SortedArray<N, [usize;2]>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize> DerefMut for ElementaryCubicPath<N>{
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}


impl<const N: usize> ElementaryCubicPath<N> {
    pub fn new() -> Self {
        Self(SortedArray::new())
    }

    pub fn act_on(&self, p: &mut [usize; N]) {
        // 'p' must be an array of (unsorted) distinct points
        debug_assert!({
            let mut q = p.clone();
            q.sort();
            q.into_iter().zip( q.into_iter().skip(1) ).all(|(a, b)| a < b )
        }, "p={p:?}");

        // in order for 'self' to be able to act on 'p', start of 'self' must be contained in 'p'
        debug_assert!(
            self.iter().all(|[s,_]| p.iter().any(|x| s==x ) ),
            "self={self:?}, p={p:?}"
        );

        // send each vertex to the goal
        for x in p {
            if let Ok(idx) = self.binary_search_by_key(x, |&[a, _]| a ) {
                *x = self[idx][1]
            };
        }
    }
    
    fn from_morse_cube(c: &MorseCube<N>) -> (Self, [[usize; N];2]) {
        assert!(c.edges.len()==1, "a cubic path can be created only if the cube is 1-dimensional, but it is {} dimensional. cube={c:?}", c.edges.len());

        let edge = c.edges[0];

        let mut start = c.vertices;
        start.add(edge[0]);

        let mut end = c.vertices;
        end.add(edge[1]);

        let out = {
            let mut out = ElementaryCubicPath::new();
            out.add_without_sort(edge);
            out
        };

        (out, [start.data, end.data])
    }


    // // This function takes two 'ElementaryCubicPath' by reference and returns a composition if the two paths are composable enclosed in Ok().
    // // Otherwise, it will return 'Err(())'.
    // // The returned path might be an identity path.
    // pub fn composed_with(&self, other: &Self) -> Result<Self, ()> {
    //     let mut out = self.0;

    //     for out_handle in &mut out {
    //         // Check whether or not the path is composable at the point.
    //         // We can use the binary search because the array is ordered by the start of the path
    //         let idx = if let Ok(idx) = other.binary_search_by_key( &out_handle[1], |[a,_]| *a ) {
    //             idx
    //         } else {
    //             return Err(())
    //         };

    //         // compute the composition
    //         *out_handle = other[idx];
    //     }

    //     Ok(Self(out))
    // }

    fn is_identity(&self) -> bool {
        self.len==0
    }
}


#[derive(Clone)]
pub struct CubicPath<const N: usize> {
    path: Vec<ElementaryCubicPath<N>>,

    // 'start' is the start position of the path. It must be sorted.
    start: [usize; N],

    // 'end' is the end position of the path. It is in general not sorted.
    end: [usize; N],
}

impl<const N: usize> Deref for CubicPath<N> {
    type Target = Vec<ElementaryCubicPath<N>>;
    fn deref(&self) -> &Self::Target {
        &self.path
    }
}

impl<const N: usize> DerefMut for CubicPath<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.path
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
    pub fn act_on(&self, p: &mut [usize; N]) {
        for e in self.iter() {
            e.act_on(p);
        }
    }

    pub fn composed_with(mut self, mut other: Self) -> Self {

        // 'self.end' must coinside with 'other.start' setwize.
        // 'self.start' is sorted, but 'self.end' is in general not sorted.
        debug_assert!(
            {
                let mut end = self.end;
                end.sort();
                end == other.start
            },
            "cannot compose two path. path1 has end = {:?}, \nbut path2 has start = {:?}.",
            self.end, other.start
        );

        self.append(&mut other);

        let mut end = self.start;
        self.act_on(&mut end);
        self.end = end;
        self.reduce_to_geodesic()
        // self
    }

    // pub fn reduce_to_geodesic(self) -> Self {
    //     // do the path reduction untill the path cannnot be reduced anymore.
    //     let mut path_length = self.len();
    //     let mut path = self.reduce();

    //     while path_length != path.len() {
    //         path_length = path.len();
    //         path = self.reduce();
    //     }
    //     path
    // }

    // pub fn reduce(self) -> Self {
    //     if self.is_empty() {return self};

    //     let mut reduced_path = vec![self[0]];
    //     reduced_path.reserve( self.len() );
        
    //     for f in self.into_iter().skip(1) {
    //         if reduced_path.is_empty() {
    //             reduced_path.push( f );
    //             continue;
    //         }

    //         // if 'reduced_path' is nonempty, then 'reduced_path' compute the composition at the end of 'reduced_path'.
    //         if let Ok(h) = reduced_path.last().unwrap().composed_with( &f ) {
    //             if h.is_identity() {
    //                 reduced_path.pop();
    //             } else {
    //                 *reduced_path.last_mut().unwrap() = h;
    //             }
    //         } else {
    //             // if 'f' cannot be composed with the last element, then simply push 'f' at the end.
    //             reduced_path.push( f );
    //         };
    //     }

    //     Self( reduced_path )
    // }

    fn reduce_to_geodesic(self) -> Self {
        // let mut out = Vec::new();

        // for mut f in self.0 {
        //     if out.is_empty() {
        //         out.push( f );
        //         continue;
        //     }

            
        //     for edge_handle in f.0.iter_mut().filter(|[x,y]| x!=y ) {
        //         let mut identity_path_idx = Vec::new();

        //         let [start, end] = *edge_handle;
                
        //         // check if the edge motion can be brought to front or not.
        //         let mut prev_handle = edge_handle;
        //         for (i, curr_path) in out.iter_mut().enumerate().rev() {
        //             if let Ok(j) = curr_path.binary_search_by_key(&start, |[s,_]| *s) {
        //                 prev_handle[0] = end;
        //                 prev_handle = &mut curr_path[j];
        //                 prev_handle[1] = end;
        //             } else {
        //                 if let Ok(j) = curr_path.binary_search_by_key(&end, |[s,_]| *s) {
        //                     // then the edge motions cancel out
        //                     prev_handle[0] = end;
        //                     curr_path[j][1] = end;

        //                     // Now 'curr_path' might be an identity. In this case remember to erase it later.
        //                     // if curr_path.is_identity() {
        //                     //     identity_path_idx.push(i);
        //                     // }
        //                 } else {
        //                     // otherwise, we leave everything as is
        //                 }
        //                 break;
        //             };
        //         }

        //         identity_path_idx.sort();
        //         debug_assert!(identity_path_idx.iter().zip(identity_path_idx.iter().skip(1)).all(|(x, y)| x!=y));
        //         for idx in identity_path_idx.into_iter().rev() {
        //             out.remove(idx);
        //         }
        //     }

        //     // Now 'f' might be changed.
        //     // add 'f' to out if 'f' is not an identity. 
        //     if !f.is_identity() {
        //         out.push(f)
        //     }
        // }
        let (start, end) = (self.start, self.end);
        let mut out = Vec::new();

        for f in self.path {
            if out.is_empty() {
                out.push(f);
                continue;
            }

            let mut new_motion = ElementaryCubicPath::new();
            for [v, w] in f.into_iter() {
                if let Some(non_commuting_idx) = (0..out.len()).rev().find(|&i| out[i].iter().any(|&[x, y]| v==x || v==y || w==x || w==y ) ) {
                    let non_commuting_motion = &mut out[non_commuting_idx];
                    if let Ok(idx) = non_commuting_motion.binary_search(&[w, v] ) {
                        // if 'non_commuting_motion' contains the inverse motion of [v,w], remove it from 'non_commuting_motion'
                        non_commuting_motion.remove(idx);
                        
                        // remove 'non_commuting_motion' if the motion is now the identity motion
                        if non_commuting_motion.is_identity() {
                            out.remove(non_commuting_idx);
                        }
                    } else {
                        // if 'non_commuting_motion' does not contain the inverse
                        if non_commuting_idx == out.len()-1 {
                            new_motion.add_without_sort([v,w]);
                        } else {
                            out[non_commuting_idx+1].add([v,w]);
                        }
                    }
                } else {
                    // if '[v,w]' commutes with every motion in 'out', then insert the motion at the front.
                    if [v, w] == [85,86] {
                        println!("out={out:?}");
                    }
                    out[0].add([v,w]);
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
}

pub struct MorsePath<'a, const N: usize> {
    // 'start' is sorted at any given time.
    start: [usize; N],

    // 'end' is in general not sorted.
    end: [usize; N],

    path: Vec<MorseCube<N>>,
    graph: &'a RawSimpleGraph,
}

impl<'a, const N: usize> MorsePath<'a, N> {
    pub fn new(start: [usize; N], end: [usize; N], graph: &'a RawSimpleGraph) -> Self {
        assert!(
            start.iter().zip(start.iter().skip(1)).all(|(x,y)| x<y ),
            "'start' must be sorted, but it is not. start = {start:?}."
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

    pub fn compose(&mut self, mut other: MorsePath<N>) {
        self.path.append(&mut other.path);
    }

    pub fn get_geodesic(&'a self) -> CubicPath<N> {
        let start_path = {
            let mut path = VecDeque::new();
            let base = get_edge_path_recc::<N, true>(&self.start, &mut path, self.graph);
            let path: Vec<_> = path.into();
            CubicPath{ path, start: self.start, end: base }
        };

        let end_path = {
            let mut path = VecDeque::new();
            // because 'self.end' is in general not sorted, we have to sort it in order to run 'get_edge_path_recc()'.
            let mut end = self.end; end.sort();
            let base = get_edge_path_recc::<N, false>(&end, &mut path, self.graph);
            let path: Vec<_> = path.into();
            CubicPath{ path, start: base, end }
        };

        println!("start_path={start_path:?}");
        println!("end_path={end_path:?}");

        self.path.iter().map(|cell| cell.get_edge_path(self.graph) )
            .chain(std::iter::once(end_path))
            .fold( start_path, |accum, x| accum.composed_with(x) )
    }
}


#[cfg(test)]
mod morse_path_test {
    // 'N' is the number of robots
    const N: usize = 5;
    use crate::{ graph_collection, graphics, MorsePath, MorseCube, SortedArray };

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


fn get_edge_path_recc<const N: usize, const FORWARD: bool>(points: &[usize; N], path: &mut VecDeque<ElementaryCubicPath<N>>, graph: &RawSimpleGraph) -> [usize; N] {
    // points must be ordered, because we will use the binary search on it.
    debug_assert!(
        points.into_iter().zip(points.into_iter().skip(1)).all(|(p, q)| p < q ),
        "points={points:?} \nhave to be ordered at this point, but it is not."
    );


    // base case
    if points.iter().enumerate().all(|(i, &x)| i==x ) {
        return *points;
    }

    let (elementary_path, next_points) = {
        let mut out = SortedArray::new();
        let mut next_points = [0; N];
        for (j, p) in points.iter().copied().enumerate() {
            let next = graph.adjacent_vertices_iter(p)
                .take_while(|&v| v < p)
                .filter( |&v| graph.maximal_tree_contains([v, p]) )
                .next().unwrap_or(p);
            
            // check whether 'next' is occupied by some point from 'points' or not.
            // if it does, then 'p' cannot proceed to 'next'.
            let next = if points.binary_search(&next).is_ok() {
                p
            } else {
                // Now 'p' can proceed to next if it is order-respecting.
                if next+1 == p {
                    next
                } else if points.iter().all(|q| !(next..p).contains(q) ) { // TODO: this can be more efficient
                    next
                } else {
                    // Otherwise 'p' cannnot proceed to 'next'.
                    p
                }
            };

            if p != next {
                if FORWARD {
                    out.add([p, next]);
                } else {
                    out.add([next, p]);
                };
            }

            // record 'next'
            next_points[j] = next;
        }
        (ElementaryCubicPath( out ), next_points)
    };

    if FORWARD {
        path.push_back( elementary_path );
    } else {
        path.push_front( elementary_path );
    }

    // println!("next_points = {next_points:?}");

    // run the function recursively
    get_edge_path_recc::<N, FORWARD>(&next_points, path, graph)
}