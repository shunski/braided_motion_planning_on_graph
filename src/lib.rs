use std::{collections::VecDeque, fmt::Debug};

use topo_spaces::graph::RawSimpleGraph;

use std::ops::{Deref, DerefMut};
use std::mem;


pub mod graph_collection;
pub mod graphics;
pub mod search;
pub mod operators;


#[derive(Clone)]
// This struct represents a cell in a discrete morse complex built upon the cube complex.
// 'N' is the number of points (or robots)
// 'D' is the dimension of the cell
pub struct MorseCube<const N: usize, const D: usize> {
    // cube is a product of DIRECTED edges
    // each edge is directed by: cube[i][0] -> cube[i][1]
    pub cube: [[usize; 2]; D],

    // number of points in the basepoint
    pub n_points_stacked_at_basepoint: usize,

    // number of points stacked on the cube.
    // 'n_points_stacked_at_cube[d][i]' is the number of points that is stacked at the 'cube[d][i]'
    pub n_points_stacked_at_cube: [[Vec<usize>; 2]; D]
}

// This strust represents the one-step of the cubic path.
// Any instance of this struct must satisfy the condition that the array 'self[..][0]' is ordered, that is, the start of elementary path is sorted. 
#[derive(Clone, Copy)]
pub struct ElementaryCubicPath<const N: usize>{
    data: [[usize;2]; N],
    size: usize,
}

impl<const N: usize> Deref for ElementaryCubicPath<N> {
    type Target = [[usize;2]];
    fn deref(&self) -> &Self::Target {
        &self.data[..self.size]
    }
}

impl<const N: usize> DerefMut for ElementaryCubicPath<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data[..self.size]
    }
}

impl<const N: usize> Debug for ElementaryCubicPath<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data[..self.size].fmt(f)
    }
}

impl<const N: usize> ElementaryCubicPath<N> {
    pub fn new() -> Self {
        Self { data: [[0; 2]; N], size: 0 }
    }

    pub fn act_on(&self, p: &mut [usize; N]) {
        println!("HELLO");
        // 'p' must be an array of (unsorted) distinct points
        debug_assert!({
            let mut q = p.clone();
            q.sort();
            q.into_iter().zip( q.into_iter().skip(1) ).all(|(a, b)| a < b )
        }, "p={p:?}");

        // in order for 'self' to be able to act on 'p', start of 'self' must be contained in 'p'
        debug_assert!(
            self.iter().all(|[s,_]| p.iter().any(|x| s==x ) )
        );

        // send each vertex to the goal
        for x in p {
            if let Ok(idx) = self.binary_search_by_key(x, |&[a, _]| a ) {
                *x = self[idx][1]
            };
        }
    }

    fn into_iter(self) -> impl Iterator<Item = [usize; 2]> { 
        self.data.into_iter().take(self.size)
    }

    fn add(&mut self, [v, w]: [usize; 2]) {
        debug_assert!(self.size < N); 
        debug_assert!(v < w);
        debug_assert!(self.iter().all(|&[x, y]| v!=x && v!=y && w!=x && w!=y ));

        let idx = self.binary_search_by_key(&v, |&[x,_]| x ).unwrap_err();
        let mut tmp = [v, w];
        for i in idx..self.size+1 {
            mem::swap(&mut tmp, &mut self[i]);
        }


        self.size += 1;
    }

    fn add_without_sort(&mut self, [v, w]: [usize; 2]) {
        debug_assert!(self.size < N); 
        debug_assert!(self.iter().all(|&[x, y]| v!=x && v!=y && w!=x && w!=y ));

        self.data[self.size] = [v,w];


        self.size += 1;
    }
    
    fn from_morse_cube(c: MorseCube<N,1>, graph: &RawSimpleGraph) -> (Self, [[usize; N];2]) {
        let edge = c.cube[0];

        let mut start = [edge[0]; N];
        let mut start_handle = start.iter_mut().skip(1); // skip the iterator by 1 because the first element is set

        // insert the points stacked at the basepoint
        for p in 0..c.n_points_stacked_at_basepoint {
            *start_handle.next().unwrap() = p;
        }
        
        // compute the positions of points stacked at the start vertex of the edge
        let start_idx = edge[0];
        let points_at_start = &c.n_points_stacked_at_cube[0][0];
        graph.adjacent_vertices_iter( start_idx )
            .skip_while(|&w| w < start_idx )
            .filter(|&w| graph.maximal_tree_contains([start_idx, w]))
            .zip( points_at_start.iter() )
            .for_each(|(w, &n_stacked_at_v)| {
                (0..n_stacked_at_v).for_each(|i| *start_handle.next().unwrap() = w+i );
            });

        // compute the positions of points stacked at the end vertex of the edge
        let end_idx = edge[1];
        let points_at_end = &c.n_points_stacked_at_cube[0][1];
        graph.adjacent_vertices_iter( end_idx )
            .skip_while(|&w| w < end_idx )
            .filter(|&w| graph.maximal_tree_contains([start_idx, w]))
            .zip( points_at_end.iter() )
            .for_each(|(w, &n_stacked_at_v)| {
                (0..n_stacked_at_v).for_each(|i| *start_handle.next().unwrap() = w+i );
            });
        
        // the two iteration above should exactly use up the iterator 'start_handle'.
        assert_eq!(start_handle.next(), None);

        // sort the elementary motion so that the start is ordered
        start.sort();

        let out = {
            let mut out = [[0;2]; N];
            out[0] = edge;
            Self{
                data: out,
                size: 1
            }
        };

        let mut end = start;
        
        out.act_on(&mut end);
        end.sort();

        (out, [start, end])
    }

    fn remove(&mut self, [v,w]: [usize; 2]) {
        debug_assert!(self.size > 0, "cannot remove from edge motion from identity motion");
        
        let idx = self.binary_search(&[v, w]).unwrap();
        for i in idx..self.size-1 {
            self.swap(i, i+1);
        }

        self.size -= 1;
    }


    fn remove_by_idx(&mut self, idx: usize) {
        debug_assert!(self.size > 0, "cannot remove from edge motion from identity motion");
        
        for i in idx..self.size-1 {
            self.swap(i, i+1);
        }

        self.size -= 1;
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
        self.size==0
    }
}


#[derive(Clone)]
pub struct CubicPath<const N: usize> (Vec<ElementaryCubicPath<N>>);

impl<const N: usize> Deref for CubicPath<N> {
    type Target = Vec<ElementaryCubicPath<N>>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize> DerefMut for CubicPath<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<const N: usize> Debug for CubicPath<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for p in self.0.iter() {
            write!(f, "{p:?}\n")?
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
        self.append(&mut other);
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

        let mut out = Vec::new();

        for f in self.0 {
            if out.is_empty() {
                out.push(f);
                continue;
            }

            let mut new_motion = ElementaryCubicPath::new();
            for [v, w] in f.into_iter() {
                if let Some(non_commuting_idx) = (0..out.len()).rev().find(|&i| out[i].iter().any(|&[x, y]| v==x || v==y || w==x || w==y ) ) {
                    let non_commuting_motion = &mut out[non_commuting_idx];
                    if let Ok(idx) = non_commuting_motion.binary_search(&[w, v] ) {
                        // if 'non_commuting_motion' contains an inverse motion of [v,w], remove it from 'non_commuting_motion'
                        non_commuting_motion.remove_by_idx(idx);
                        
                        // remove 'non_commuting_motion' if the motion is now the identity motion
                        if non_commuting_motion.is_identity() {
                            out.remove(non_commuting_idx);
                        }
                    } else {
                        if non_commuting_idx == out.len()-1 {
                            new_motion.add_without_sort([v,w]);
                        } else {
                            out[non_commuting_idx+1].add([v,w])
                        }
                    }
                } else {
                    // if '[v,w]' commutes with every motion in 'out', then insert the motion at the front.
                    out[0].add([v,w]);
                }
            }

            // add 'new_motion' to 'out' if 'new_motion' is not the identity motion
            if !new_motion.is_identity() {
                out.push(new_motion);
            }

        }

        Self( out )
    }
}

impl<const N: usize, const D: usize> MorseCube<N, D> {
    pub fn is_valid(&self, graph: &RawSimpleGraph) -> bool {
        // check that the number of points is correct
        let n_stacked_points = self.n_points_stacked_at_cube.iter().map(|[x1,x2]| x1.iter().chain(x2.iter()) ).flatten().sum::<usize>();
        if N != D+self.n_points_stacked_at_basepoint+n_stacked_points {
            return false;
        }

        // check that the cell is contained in the graph
        self.cube.iter()
            .zip(self.n_points_stacked_at_cube.iter())
            .map(|([x1,x2], [y1,y2])| [(x1, y1), (x2, y2)] )
            .flatten()
            .map(|(v, x)| graph.adjacent_vertices_iter(*v).skip_while(|&w| w < *v ).filter(|&w| graph.maximal_tree_contains([*v, w])).zip( x.iter() ) )
            .flatten()
            .all(|(next_vertex, &n_points) | 
                // check that none of them is neither a vertex of degree one nor an essential vertex
                (0..n_points-1).all(|i| graph.degree_of(next_vertex + i)==2 )
            )
    }
}

impl<const N: usize> MorseCube<N, 1> {
    pub fn get_edge_path(self, graph: &RawSimpleGraph) -> CubicPath<N> {
        let (critical_motion, [start, end]) = ElementaryCubicPath::from_morse_cube(self, graph);

        // 'start' and 'end' are already sorted.
        debug_assert!( 
            end.into_iter().zip(end.into_iter().skip(1)).all(|(a,b)| a < b ),
            "end = {end:?} \n is not ordered. implimentation of 'ElementaryCubicPath::from_morse_cube' might be wrong."
        );
        debug_assert!( 
            start.into_iter().zip(start.into_iter().skip(1)).all(|(a,b)| a < b ),
            "start = {start:?} \n is not ordered. implimentation of 'ElementaryCubicPath::from_morse_cube' might be wrong."
        );

        let mut path = VecDeque::from([critical_motion]);
        // panic!("start={start:?}");
        // panic!("end={end:?}");
        Self::get_edge_path_recc::<false>(start, &mut path, graph);
        Self::get_edge_path_recc::<true>(end, &mut path, graph);

        let path = CubicPath(
            path.into()
        );

        path
    }

    fn get_edge_path_recc<const FORWARD: bool>(points: [usize; N], path: &mut VecDeque<ElementaryCubicPath<N>>, graph: &RawSimpleGraph) {
        // points must be ordered, because we will use the binary search on it.
        debug_assert!(
            points.into_iter().zip(points.into_iter().skip(1)).all(|(p, q)| p < q ),
            "points={points:?} \nhave to be ordered at this point, but it is not."
        );


        // base case
        if points.iter().enumerate().all(|(i, &x)| i==x ) {
            return;
        }

        let (cubic_path, size, next_points) = {
            let mut out = [[0;2]; N];
            let mut next_points = [0; N];
            let mut i = 0;
            for (j, p) in points.into_iter().enumerate() {
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
                    out[i] = if FORWARD {
                        [p, next]
                    } else {
                        [next, p]
                    };
                    i+=1;
                }

                // record 'next'
                next_points[j] = next;
            }
            (out, i, next_points)
        };

        let elementary_path = ElementaryCubicPath{ data: cubic_path, size };

        if FORWARD {
            path.push_back( elementary_path );
        } else {
            path.push_front( elementary_path );
        }

        // println!("next_points = {next_points:?}");

        // run the function recursively
        Self::get_edge_path_recc::<FORWARD>(next_points, path, graph)
    }
}

// fn get_a_path_on_maximal_tree<const N: usize>(start: [usize;N], goal: [usize;N], graph: &RawSimpleGraph) -> Vec<MorseCube<N,1>> {
//     s
// }

fn get_a_path_on_maximal_tree_recc<const N: usize>(start: &[usize;N], goal: &[usize;N], essential_vertex: usize, graph: &RawSimpleGraph) -> Vec<MorseCube<N,1>> {
    let out = Vec::new();
    
    let next_vertices = graph.adjacent_vertices_iter(essential_vertex)
        .skip_while(|v| *v < essential_vertex)
        .filter(|v| graph.maximal_tree_contains([essential_vertex, *v]) )
        .collect::<Vec<_>>();

    let mut next_essential_vertices = Vec::new(); 
    next_essential_vertices.reserve( next_vertices.len() );
    for v in &next_vertices {
        let next_essentail_vertex = (*v..).find(|w| graph.degree_of(*w) > 2 );
        next_essential_vertices.push( next_essentail_vertex );
    }

    // collect points below the essential vertex
    let points_start_here = (0..=essential_vertex).rev().filter(|w| graph.degree_of(*w)==2 ).filter(|w| start.binary_search(w).is_ok() ).collect::<Vec<_>>();
    let points_goal_here = (0..=essential_vertex).rev().filter(|w| graph.degree_of(*w)==2 ).filter(|w| goal.binary_search(w).is_ok() ).collect::<Vec<_>>();

    

    out
}


// #[cfg(test)]
// mod cubic_path_test {
//     use crate::ElementaryCubicPath;

//     #[test] 
//     fn reduce_to_geodesic() {
//         let a = ElementaryCubicPath([[0,0], [1,2]]);
//     }
// }