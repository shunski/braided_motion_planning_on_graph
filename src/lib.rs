use std::collections::VecDeque;

use topo_spaces::graph::RawSimpleGraph;

use std::ops::{Deref, DerefMut};



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
// Any instance of this struct must satisfy the condition that the array 'self.0[..][0]' is ordered, that is, the start of elementary path is sorted. 
#[derive(Clone, Copy)]
pub struct ElementaryCubicPath<const N: usize>([[usize;2]; N]);

impl<const N: usize> Deref for ElementaryCubicPath<N> {
    type Target = [[usize;2]; N];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<const N: usize> DerefMut for ElementaryCubicPath<N> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}


impl<const N: usize> ElementaryCubicPath<N> {
    pub fn between(start: [usize; N], end: [usize; N]) -> Self {
        let mut critical_motion = [[0;2];N];
        for i in 0..N {
            critical_motion[i] = [start[i], end[i]];
        }
        Self(critical_motion)
    }

    pub fn act_on(&self, p: &mut [usize; N]) {
        // 'p' is the set of distinct vertices on the graph.
        // though 'p' is not ordered in general, 'p' has to coinside start of 'self' setwise.
        debug_assert!({
            let mut q = p.clone();
            q.sort();
            q.into_iter().zip( self.into_iter() ).all(|(a,[b,_])| a==b )
        });

        // send each vertex to the goal
        for x in p {
            let idx = self.binary_search_by_key(x, |&[a, _]| a ).unwrap();
            *x = self[idx][1];
        }
    }
    
    fn from_morse_cube(c: MorseCube<N,1>, graph: &RawSimpleGraph) -> Self {
        let edge = c.cube[0];

        let mut out = [edge; N];
        let mut out_handle = out.iter_mut().skip(1); // skip the iterator by 1 because the first element is set

        // insert the points stacked at the basepoint
        for p in 0..c.n_points_stacked_at_basepoint {
            *out_handle.next().unwrap() = [p; 2];
        }
        
        // compute the positions of points stacked at the start vertex of the edge
        let start_idx = edge[0];
        let points_at_start = &c.n_points_stacked_at_cube[0][0];
        graph.adjacent_vertices_iter( start_idx )
            .skip_while(|&w| w < start_idx )
            .filter(|&w| graph.maximal_tree_contains([start_idx, w]))
            .zip( points_at_start.iter() )
            .for_each(|(w, &n_stacked_at_v)| {
                (0..n_stacked_at_v).for_each(|i| *out_handle.next().unwrap() = [w+i;2] );
            });

        // compute the positions of points stacked at the end vertex of the edge
        let end_idx = edge[1];
        let points_at_end = &c.n_points_stacked_at_cube[0][1];
        graph.adjacent_vertices_iter( end_idx )
            .skip_while(|&w| w < end_idx )
            .filter(|&w| graph.maximal_tree_contains([start_idx, w]))
            .zip( points_at_end.iter() )
            .for_each(|(w, &n_stacked_at_v)| {
                (0..n_stacked_at_v).for_each(|i| *out_handle.next().unwrap() = [w+i;2] );
            });
        
        // the two iteration above should exactly use up the iterator 'out_handle'.
        assert_eq!(out_handle.next(), None);

        // sort the elementary motion so that the start is ordered
        out.sort_by_key(|k| k[0] );

        Self(out)
    }


    // This function takes two 'ElementaryCubicPath' by reference and returns a composition if the two paths are composable enclosed in Ok().
    // Otherwise, it will return 'Err(())'.
    // The returned path might be an identity path.
    pub fn composed_with(&self, other: &Self) -> Result<Self, ()> {
        let mut out = self.0;

        for out_handle in &mut out {
            // Check whether or not the path is composable at the point.
            // We can use the binary search because the array is ordered by the start of the path
            let idx = if let Ok(idx) = other.binary_search_by_key( &out_handle[1], |[a,_]| *a ) {
                idx
            } else {
                return Err(())
            };

            // compute the composition
            *out_handle = other[idx];
        }

        Ok(Self(out))
    }

    fn is_identity(&self) -> bool {
        self.iter().all(|[x,y]| x==y )
    }
}

impl<const N: usize> std::convert::From<[[usize;2]; N]> for ElementaryCubicPath<N> {
    fn from(value: [[usize;2]; N]) -> Self {
        Self(value)
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

impl<const N: usize> CubicPath<N> {
    pub fn act_on(&self, p: &mut [usize; N]) {
        for e in self.iter() {
            e.act_on(p);
        }
    }

    pub fn composed_with(mut self, mut other: Self) -> Self {
        self.append(&mut other);
        self.reduce_to_geodesic()
    }

    pub fn reduce_to_geodesic(self) -> Self {
        // do the path reduction untill the path cannnot be reduced anymore.
        let mut path_length = self.len();
        let mut path = self.reduce();

        while path_length != path.len() {
            path_length = path.len();
            path = self.reduce();
        }
        path
    }

    pub fn reduce(self) -> Self {
        if self.is_empty() {return self};

        let mut reduced_path = vec![self[0]];
        reduced_path.reserve( self.len() );
        
        for f in self.into_iter().skip(1) {
            if reduced_path.is_empty() {
                reduced_path.push( f );
                continue;
            }

            // if 'reduced_path' is nonempty, then 'reduced_path' compute the composition at the end of 'reduced_path'.
            if let Ok(h) = reduced_path.last().unwrap().composed_with( &f ) {
                if h.is_identity() {
                    reduced_path.pop();
                } else {
                    *reduced_path.last_mut().unwrap() = h;
                }
            } else {
                // if 'f' cannot be composed with the last element, then simply push 'f' at the end.
                reduced_path.push( f );
            };
        }

        Self( reduced_path )
    }

    fn reduce_to_geodesic_(self) -> Self {
        let mut out = Vec::new();

        for f in self.0 {
            if out.is_empty() {
                out.push( f );
                continue;
            }

            for edge_handle in f.0.iter_mut().filter(|[x,y]| x!=y ) {
                // check if the edge motion can be brought to front or not.
                
            }
        }

        out
    }

    pub fn inverse() -> Self {
        s
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
        let critical_motion = ElementaryCubicPath::from_morse_cube(self, graph);

        let mut start = [0; N];
        let mut end = [0; N];
        let pos_iter_mut = start.iter_mut().zip(end.iter_mut());
        for ([x, y], (start_mut, end_mut)) in critical_motion.0.into_iter().zip( pos_iter_mut ) {
            *start_mut = x;
            *end_mut = y;
        }

        // order the 'end'. Note that 'start' is already ordered
        end.sort();
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

        let (cubic_path, next_points) = {
            let mut out = [[0;2]; N];
            let mut next_points = [0; N];
            for (i, &p) in points.iter().enumerate() {
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

                out[i] = if FORWARD {
                    [p, next]
                } else {
                    [next, p]
                };

                // record 'next'
                next_points[i] = next;
            }
            (out, next_points)
        };

        let elementary_path = ElementaryCubicPath::from(cubic_path);

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

fn get_a_path_on_maximal_tree<const N: usize>(start: [usize;N], goal: [usize;N], graph: &RawSimpleGraph) -> Vec<MorseCube<N,1>> {
    s
}

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

