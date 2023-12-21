use std::cmp::Ordering;

use quasi_iter::merge_sorted::merge_sorted_dyn_by_key;
use topo_spaces::graph::RawSimpleGraph;

use crate::{MorsePath, MorseCube, SortedArray};

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


pub fn path_in_tree<'a, const N: usize>(points: &mut [[usize; 2]; N], graph: &'a RawSimpleGraph) -> MorsePath<'a, N> {
    let graph = GraphInformation::<'a>::from(graph);

    let first_essential_vertex = *graph.next_essential_vertices[0].first().unwrap();
    path_at_essential_vertex_dyn(first_essential_vertex, points, &graph)
}

fn path_at_essential_vertex_dyn<'a, 'b, const N: usize>(essential_vertex: usize, points: &mut [[usize; 2]; N], graph: &'b GraphInformation<'a>) -> MorsePath<'a, N> 
{
    let (start, end): (Vec<_>, Vec<_>) = points.iter().copied().map(|[x,y]| (x,y) ).unzip();
    let start = start.try_into().unwrap();
    let end = end.try_into().unwrap();

    // initialize 'motions', which is the output of this function.
    let mut motions = MorsePath::<'a, N>::new(start, end, graph.graph_ref);
    
    motions.compose( sort_points_in_stem::<N>(essential_vertex, points, graph) );

    for &next_essential_vertex in &graph.next_essential_vertices[essential_vertex] {
        motions.compose( path_at_essential_vertex_dyn(next_essential_vertex, points, graph) );
    }

    // bring all the points to the basepoint
    let ordered_iter = collect_sorted_points(*graph.essential_vertices.first().unwrap(), points, graph);
    for (n_stacked_at_basepoint, [mut x, _]) in ordered_iter.enumerate() {
        let vertices = points.iter().map(|[s,_]| *s ).collect::<Vec<_>>().try_into().unwrap();
        while x != n_stacked_at_basepoint {
            // let 'x' flow by one step.
            let y = graph.adjacent_vertices_iter(x).filter(|&y| graph.maximal_tree_contains([x, y]) ).next().unwrap();
            let edge = SortedArray::try_from([[x,y]]).unwrap();

            // add the 1-dimensional cube if it is a critical cell.
            if let Some(morse_cube) = MorseCube::new_checked(edge, vertices) {
                motions.add( morse_cube );
            }

            // update 'x'
            x = y;
        }
    }

    motions
}

fn sort_points_in_stem<'a, const N: usize>(essential_vertex: usize, points: &mut [[usize;2]; N], graph: &'a GraphInformation) -> MorsePath<'a, N> {
    let (start, end): (Vec<_>, Vec<_>) = points.iter().copied().map(|[x,y]| (x,y) ).unzip();
    let start = start.try_into().unwrap();
    let end = end.try_into().unwrap();
    let mut motions: MorsePath<'_, N> = MorsePath::new(start, end, graph);

    // 'next_branch' is the smallest vertex that is greater than any child vertices of the essential vertex.
    let next_branch = {
        let prev_essential_vertex = (0..essential_vertex).rev()
            .find(|&idx| graph.degree_in_maximal_tree(idx) != 2 )
            .unwrap();
        if prev_essential_vertex == 0 {
            0
        } else {
            *graph.next_vertices[prev_essential_vertex].iter().find(|&v| v > &essential_vertex).unwrap_or(&usize::MAX)
        }
    };
    
    let mut points_in_stem_idx = (0..essential_vertex).rev()
        .take_while(|&i| graph.degree_in_maximal_tree(i)!=2 )
        .filter(|i| points.binary_search_by_key(i, |p| p[0] ).is_ok() )
        .collect::<Vec<_>>();


    let is_travelling = |idx: usize, points: &[[usize;2];N]| -> bool {
        points[idx][1] <= essential_vertex || points[idx][1] >= next_branch
    };

    let cmp = |idx1: usize, idx2: usize, points: &[[usize;2];N]| -> Ordering {
        if is_travelling(idx1, points) && !is_travelling(idx2, points) {
            Ordering::Less
        } else if !is_travelling(idx1, points) && is_travelling(idx2, points) {
            Ordering::Greater
        } else {
            // otherwise the ordering is determined by the goal vertex
            points[idx1][1].cmp( &points[idx2][1] )
        }
    };

    let needs_swapping = |idx: usize, points: &[[usize;2];N]| -> bool {
        (0..idx).filter(|&i| !is_travelling(i, points))
            .filter(|&i| points[i][1] > points[idx][1] )
            .any(|i| 
                // the folloing expression returns true if the goal of 'i' and the goal 'idx' belongs to the same branch
                graph.next_vertices[i]
                    .iter()
                    .zip( graph.next_vertices[i].iter().skip(1).chain(&[next_branch])) 
                    .any( |(&i, &j)| (i..j).contains(&points[i][1]) && (i..j).contains(&points[idx][1]) )
            )
    };
    

    // The base case: if 'points_in_stem' are already sorted, then return.
    if points_in_stem_idx.iter()
        .zip(points_in_stem_idx.iter().skip(1))
        .all(|(&i, &j)| cmp(i, j, points) == Ordering::Less ) 
    {
        return motions;
    }


    while !points_in_stem_idx.is_empty() && !needs_swapping( *points_in_stem_idx.last().unwrap(), points) {
        let idx = points_in_stem_idx.pop().unwrap();
        let branch = *graph.next_vertices[essential_vertex].iter()
            .take_while(|&&i| i < points[idx][1] )
            .last().unwrap();

        // compute the motion and record it.
        motions.add( push::<true, N>([essential_vertex, branch], idx, points) );
    }

    // Having processed the trivial moves, do the swapping at the essential vertex (only once)
    let n_branches = graph.next_vertices[essential_vertex].len();
    let chunks = {
        let mut cuts = points_in_stem_idx.iter().rev().skip(1)
            .zip( points_in_stem_idx.iter().rev() )
            .filter(|(&x, &y)| cmp(x, y, points) == Ordering::Greater )
            .map(|(&x,_)| x )
            .take(n_branches)
            .collect::<Vec<_>>();
        cuts.reverse();
        cuts.push(usize::MAX);
        cuts.iter().zip(cuts.iter().skip(1)).map(|(&x, &y)| x..y ).collect::<Vec<_>>()
    };
    debug_assert!(chunks.len() <= n_branches);

    // push all the points to above the essential vertex
    for (chunk, &next_vertex) in chunks.iter().rev().zip(graph.next_vertices[essential_vertex].iter().rev()) {
        for i in chunk.clone().into_iter().rev() {
            push::<true, N>([essential_vertex, next_vertex], points[i][0], points);
        }
    }

    let mut falling_points = chunks.iter()
        .map(|r| r.len() )
        .enumerate()
        .map(|(i, l)| {
            let next_vertex = graph.next_vertices[essential_vertex][i];
            (next_vertex..next_vertex+l, next_vertex)
        })
        .map(|(mut r, next_vertex)| (r.next().unwrap(), r, next_vertex) )
        .collect::<Vec<_>>();

    falling_points.sort_by(|(v,_,_), (w,_,_)| {
        let idx1 = points.binary_search_by_key(v, |[s,_]| *s ).unwrap();
        let idx2 = points.binary_search_by_key(w, |[s,_]| *s ).unwrap();
        cmp(idx1, idx2, points) 
    });

    while !falling_points.is_empty() {
        let (v, iter, next_vertex) = falling_points.first_mut().unwrap();
        let idx = points.binary_search_by_key(v, |[s,_]| *s ).unwrap();
        motions.add( push::<false, N>([*next_vertex, essential_vertex], idx, points));
        if let Some(w) = iter.next() {
            *v = w;

            // sort the 'falling_points' again.
            // ToDo: This sort can be more efficient.
            falling_points.sort_by(|(v,_,_), (w,_,_)| {
                let idx1 = points.binary_search_by_key(v, |[s,_]| *s ).unwrap();
                let idx2 = points.binary_search_by_key(w, |[s,_]| *s ).unwrap();
                cmp(idx1, idx2, points) 
            });
        } else {
            // if 'iter.next()' return 'None', then remove the tuple from 'falling_points'.
            falling_points.remove(0);
        }
    }

    motions.compose( sort_points_in_stem(essential_vertex, points, graph) );

    motions
}


// This function takes a point by 'idx: usize' and its motion by '[usize; 2]' and computes the motion'.
fn push<const UPWARD: bool, const N: usize>( edge: [usize; 2], idx: usize, points: &mut [[usize; 2]; N] ) -> MorseCube<N> {
    let terminal = edge[1];

    // if 'terminal' is occupied by another point, then push them, and make some space.
    if let Ok(i) = points.binary_search_by_key(&terminal, |[s,_]| *s ) {
        let last_idx = if UPWARD {
            (terminal+1..).zip(i+1..).find(|&(v,j)| v != points[j][0] ).unwrap().1
        } else {
            (0..terminal).rev().zip((0..i).rev()).find(|&(v,j)| v != points[j][0] ).unwrap().1
        };

        if UPWARD {
            for j in i..last_idx { points[j][0] += 1; } 
        } else {
            for j in last_idx+1..=i { points[j][0] -= 1; }
        }
    }


    // create the motion
    let vertices: SortedArray<N, usize> = points.iter()
        .map(|&[s,_]| s )
        .filter(|&s| s!=points[idx][0])
        .collect::<Vec<_>>()
        .try_into().unwrap();
    debug_assert!(vertices.len() == N-1);

    let edges: SortedArray<N, [usize; 2]> = [edge].into();

    let motion = MorseCube::new_unchecked(edges, vertices);


    // finally, physically move the point
    points[idx][0] = terminal;
    points.sort(); // ToDo: This can be more efficient

    motion
}

fn collect_sorted_points<const N: usize>(essential_vertex: usize, points: &[[usize; 2]; N], graph: &GraphInformation) -> Box<dyn Iterator<Item = [usize; 2]>> {
    let points_in_stem = (0..essential_vertex).rev()
        .take_while(|v| v==&essential_vertex || graph.degree_in_maximal_tree(*v)==2 )
        .filter_map(|v| points.binary_search_by_key(&v, |[s,_]| *s ).ok())
        .map(|idx| points[idx] )
        .collect::<Vec<_>>();

    let next_iters = graph.next_essential_vertices[essential_vertex].iter()
        .map( |&v| 
            collect_sorted_points(v, points, graph)
        );

    let next_iters_merged = merge_sorted_dyn_by_key( next_iters, |[s,_]| *s );

    Box::new(
        points_in_stem.into_iter().chain( next_iters_merged )
    )

}

struct UsizeIndexable<T>(Vec<(usize, T)>);

impl<T> std::ops::Index<usize> for UsizeIndexable<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[self.0.binary_search_by_key(&index, |&(i,_)| i ).unwrap()].1
    }
}

struct GraphInformation<'a> {
    pub graph_ref: &'a RawSimpleGraph,
    pub next_essential_vertices: UsizeIndexable<Vec<usize>>,
    pub next_vertices: UsizeIndexable<Vec<usize>>,
    pub essential_vertices: Vec<usize>, // vertices of degree not equal to 2 in maximal tree
}

impl<'a> std::ops::Deref for GraphInformation<'a> {
    type Target = RawSimpleGraph;
    fn deref(&self) -> &Self::Target {
        &self.graph_ref
    }
}

impl<'a> GraphInformation<'a> {
    pub fn from(graph_ref: &'a RawSimpleGraph) -> Self {
        Self {
            graph_ref,
            next_essential_vertices: Self::get_next_essential_vertices_dictionary(graph_ref),
            next_vertices: Self::get_next_vertices_dictionary(graph_ref),
            essential_vertices: Self::get_essential_vertices(graph_ref),
        }
    }

    fn get_next_essential_vertices_dictionary(graph: &RawSimpleGraph) -> UsizeIndexable<Vec<usize>> {
        let out = graph.vertex_iter()
            .map(|v| v.vertex() )
            .filter(|v| graph.degree_in_maximal_tree(*v) != 2 )
            .map(|v| (v, get_next_essential_vertices(v, graph)))
            .collect::<Vec<_>>();

        UsizeIndexable(out)
    }
    
    
    fn get_next_vertices_dictionary(graph: &RawSimpleGraph) -> UsizeIndexable<Vec<usize>> {
        let out = graph.vertex_iter()
            .map(|v| v.vertex() )
            .filter(|v| graph.degree_in_maximal_tree(*v) != 2 )
            .map(|v| (v, get_next_vertices(v, graph)))
            .collect::<Vec<_>>();

        UsizeIndexable(out)
    }
    
    fn get_essential_vertices(graph: &RawSimpleGraph) -> Vec<usize> {
        graph.vertex_iter()
        .map(|v| v.vertex() )
        .filter(|v| graph.degree_in_maximal_tree(*v) != 2 )
        .collect::<Vec<_>>()
    }
}


#[cfg(test)]
mod path_generation_on_tree_tests {
    // 'N' is the number of robots
    const N: usize = 5;
    use crate::{ graph_collection, graphics, search::GraphInformation };

    #[test]
    fn sort_points_in_stem() -> Result<(), Box<dyn std::error::Error>> {
        let (graph, embedding, name) = graph_collection::RawGraphCollection::Grid.get();
        let graph = GraphInformation::from(&graph);

        let mut points = [[0,5], [16,20], [17,19], [21,21], [49,49]];

        let morse_path = super::sort_points_in_stem(20, &mut points, &graph);

        let path = morse_path.get_geodesic();

        // draw the generated paths
        println!("drawing paths...");
        graphics::draw_edge_path::<N>(&path, &"sort_points_in_stem_test", &name, &embedding, &graph)?;

        Ok(())
    }
}

