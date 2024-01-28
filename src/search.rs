use std::cmp::Ordering;

use quasi_iter::merge_sorted::{merge_sorted_dyn_by_key, merge_sorted_dyn_by, MergeSortedDynBy};

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



// The main algorithm that computes the path 
pub fn path_in_tree<'a, const N: usize>(points: &mut [[usize; 2]; N], graph: &'a RawSimpleGraph) -> MorsePath<'a, N> {
    let graph = GraphInformation::<'a>::from(graph);


    // the dynamic programming and recursively sort the points in the graph at each essential vertices.
    let first_essential_vertex = *graph.next_essential_vertices[0].first().unwrap();
    let mut motions = path_at_essential_vertex_dyn(first_essential_vertex, points, &graph);

    // bring all the points to the basepoint
    let sorted_iter = collect_sorted_points(*graph.essential_vertices.first().unwrap(), points, &graph);
    for (n_stacked_at_basepoint, [mut x, _]) in sorted_iter.enumerate() {
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


// 'fn path_at_essential_vertex_dyn' implements the definition of the dynamic programming.
fn path_at_essential_vertex_dyn<'a, 'b, const N: usize>(
    essential_vertex: usize, points: &mut [[usize; 2]; N], graph: &'b GraphInformation<'a>
) -> MorsePath<'a, N>
{
    // initialize 'motions', which is the output of this function.
    let (start, end): (Vec<_>, Vec<_>) = points.iter().copied().map(|[x,y]| (x,y) ).unzip();
    let start = start.try_into().unwrap();
    let end = end.try_into().unwrap();
    let mut motions = MorsePath::<'a, N>::new(start, end, graph.graph_ref);

    // sort the points in stem so that all travelling points stay in the stem (but not necessarily sorted)
    // and that all the returning points are pushed to the child branches.
    motions.compose( sort_points_in_stem::<N, true>(essential_vertex, points, graph) );

    // Sort the travelling points in the stem and send those travelling points that need to be sorted to the child branches.
    motions.compose( sort_travelling_points::<N, true>(essential_vertex, points, graph) );
    
    
    // recursively run the the algorithm in the child branches.
    for &next_essential_vertex in &graph.next_essential_vertices[essential_vertex] {
        motions.compose( path_at_essential_vertex_dyn(next_essential_vertex, points, graph) );
    }

    motions
}

fn sort_points_in_stem<'a, const N: usize, const FORWARD: bool>(essential_vertex: usize, points: &mut [[usize;2]; N], graph: &'a GraphInformation) -> MorsePath<'a, N> {
    let (start, end): (Vec<_>, Vec<_>) = points.iter().copied().map(|[x,y]| (x,y) ).unzip();
    let start = start.try_into().unwrap();
    let end = end.try_into().unwrap();
    let mut motions: MorsePath<'_, N> = MorsePath::new(start, end, graph);

    // 'next_branch' is the smallest vertex that is greater than any child vertices of the essential vertex.
    let next_branch = {
        let prev_essential_vertex = (0..essential_vertex).rev()
            .find(|&idx| graph.degree_in_maximal_tree(idx) != 2 )
            .unwrap();
        
        *graph.next_vertices[prev_essential_vertex].iter().find(|&v| v > &essential_vertex).unwrap_or(&usize::MAX)
    };

    let is_travelling = |idx: usize, points: &[[usize;2]; N]| -> bool {
        points[idx][if FORWARD{1} else {0}] <= essential_vertex || points[idx][if FORWARD {1} else {0}] >= next_branch
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


    // this closure returns 'true' if the swapping has to be done at the current essential vertex.
    // That is, if this closure returns 'false', then we can simply push this vertex to the next child branch.
    let needs_swapping = |idx: usize, points: &[[usize;2];N]| -> bool {
        (0..idx).filter(|&i| cmp(i, idx, points).is_gt())
            .any(|i|
                // the following expression returns ture if 
                // (1) 'points[idx]' is travelling (this automatically implies 'points[i]' is also travelling), or if
                is_travelling(idx, points)
                ||
                // (2) the following expression returns true if the goal of 'i' and the goal 'idx' belong to the same child branch
                graph.next_vertices[essential_vertex]
                    .iter()
                    .zip( graph.next_vertices[essential_vertex].iter().skip(1).chain(&[next_branch])) 
                    .any( |(&a, &b)| (a..b).contains(&points[i][1]) && (a..b).contains(&points[idx][1]) )
            )
    };


    // 'points_in_stem_idx' is the points located at vertices smaller than or equal to the current essential vertex and 
    // larger than the previous essential vertex of valency greater than 3 or greater than or equal to the previous
    // essential vertex of valency 1.
    let mut points_in_stem_idx = (0..=essential_vertex).rev()
        .take_while(|&i| i==essential_vertex || graph.degree_in_maximal_tree(i) <= 2 )
        .filter_map(|i| points.binary_search_by_key(&i, |p| p[0] ).ok() )
        .collect::<Vec<_>>();
    points_in_stem_idx.reverse();
    

    println!("points_in_stem_idx={points_in_stem_idx:?}");

    while !points_in_stem_idx.is_empty() && !needs_swapping( *points_in_stem_idx.last().unwrap(), points) && !is_travelling( *points_in_stem_idx.last().unwrap(), points) {
        let idx = points_in_stem_idx.pop().unwrap();
        let branch = *graph.next_vertices[essential_vertex].iter()
            .take_while(|&&i| i <= points[idx][1] )
            .last().unwrap();

        // compute the motion and record it.
        motions.add( push::<true, N>([essential_vertex, branch], idx, points) );
    }


    // The base case: if 'points_in_stem' is already sorted and if all such points are travelling, then return.
    if points_in_stem_idx.iter()
        .zip(points_in_stem_idx.iter().skip(1))
        .all(|(&i, &j)| cmp(i, j, points).is_lt() ) 
        &&
        points_in_stem_idx.iter().all(|&i| is_travelling(i, points) )
    {
        return motions;
    }

    // Having processed the trivial moves, do the swapping at the essential vertex (only once)
    let n_branches = graph.next_vertices[essential_vertex].len();
    let swapping_chunks = {
        let mut cuts = points_in_stem_idx.iter().rev().skip(1)
            .zip( points_in_stem_idx.iter().rev() )
            .filter(|(&x, &y)| cmp(x, y, points).is_gt() )
            .map(|(_,&y)| y )
            .take(n_branches)
            .collect::<Vec<_>>();
        
        if cuts.len() < n_branches {
            cuts.push(0);
        }
        cuts.reverse();
        cuts.push( points_in_stem_idx.len() );
        cuts.iter().zip(cuts.iter().skip(1)).map(|(&x, &y)| x..y ).collect::<Vec<_>>()
    };
    debug_assert!( 1 < swapping_chunks.len() && swapping_chunks.len() <= n_branches, "swapping_chunks={swapping_chunks:?}" );
    debug_assert!(
        swapping_chunks.clone().into_iter().all(|chunk| chunk.end == points_in_stem_idx.len() || cmp(chunk.end-1, chunk.end, points).is_gt() ),
        "swapping_chunks={swapping_chunks:?}, and points are {:?}",
        points_in_stem_idx.iter().map(|&idx| points[idx] ).collect::<Vec<_>>()
    );
    debug_assert!( !(swapping_chunks.len() < n_branches) || swapping_chunks[0].end == 0 );

    // push all the points to above the essential vertex
    for (chunk, &next_vertex) in swapping_chunks.iter().rev().zip(graph.next_vertices[essential_vertex][..swapping_chunks.len()].iter().rev()) {
        for i in chunk.clone().into_iter().rev() {
            push::<true, N>([essential_vertex, next_vertex], i, points);
        }
    }

    let curr_points = points.clone();

    let falling_points: MergeSortedDynBy<usize, _> = merge_sorted_dyn_by(
        swapping_chunks.iter()
            .map(|r| r.len() )
            .enumerate()
            .map(|(i, l)| {
                let next_vertex = graph.next_vertices[essential_vertex][i];
                let out: Box::<dyn Iterator<Item=usize>> = Box::new(next_vertex..next_vertex+l);
                out
            }),
            |v, w| {
                let idx1 = curr_points.binary_search_by_key(v, |[s,_]| *s ).unwrap();
                let idx2 = curr_points.binary_search_by_key(w, |[s,_]| *s ).unwrap();
                cmp(idx1, idx2, &curr_points)
            }
    );

    for v in falling_points {
        let edge_initial = *graph.next_vertices[essential_vertex].iter().take_while(|&&w| w <= v ).last().unwrap(); // ToDo: use binary search here
        
        let idx = points.binary_search_by_key(&v, |[s,_]| *s ).unwrap();
        motions.add( push::<false, N>([edge_initial, essential_vertex], idx, points));
    }

    motions.compose( sort_points_in_stem::<N, FORWARD>(essential_vertex, points, graph) );

    motions
}


// This function takes a point by 'idx: usize' and its motion by 'edge: [usize; 2]' and computes the "push" motion.
fn push<const UPWARD: bool, const N: usize>( edge: [usize; 2], idx: usize, points: &mut [[usize; 2]; N] ) -> MorseCube<N> {
    println!("pushing along {edge:?}");
    let terminal = edge[1];

    // if 'terminal' is occupied by another point, then push them, and make some space.
    if let Ok(i) = points.binary_search_by_key(&terminal, |[s,_]| *s ) {
        let last_idx = if UPWARD {
            (terminal+1..).zip(i+1..).find(|&(v,j)| j==points.len() || v != points[j][0] ).unwrap().1
        } else {
            (0..=terminal).rev().zip((0..=i).rev()).take_while(|&(v,j)| v == points[j][0] ).last().unwrap().1
        };

        if UPWARD {
            for j in i..last_idx { points[j][0] += 1; } 
        } else {
            for j in last_idx..=i { points[j][0] -= 1; }
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


// 'fn sort_travelling_points' sorts points in the stem, and send those points that need to be sorted
// in the child branches to the child branches.
// This function requires that each point in 'points' inside the stem needs to be a travelling point. 
fn sort_travelling_points<'a, const N: usize, const FORWARD: bool>( essential_vertex: usize, points: &mut [[usize; 2]; N], graph: &'a GraphInformation ) -> MorsePath<'a, N> {
    // initialize 'motions', which is the return value of this function.
    let start = points.iter().map(|&[s,_]| s).collect::<Vec<_>>().try_into().unwrap();
    let end = points.iter().map(|&[_,g]| g).collect::<Vec<_>>().try_into().unwrap();
    let mut motions = MorsePath::new(start, end, &graph);
    
    // 'next_branch' is the smallest vertex that is greater than any child vertices of the essential vertex.
    let next_branch = {
        let prev_essential_vertex = (0..essential_vertex).rev()
            .find(|&idx| graph.degree_in_maximal_tree(idx) != 2 )
            .unwrap();
        
        *graph.next_vertices[prev_essential_vertex].iter().find(|&v| v > &essential_vertex).unwrap_or(&usize::MAX)
    };

    // find the point in the subtree that has the .
    let smallest_pt_idx = {
        
    };
    
    // collect the points in stem.
    let points_in_stem_idx = (0..=essential_vertex).rev()
        .take_while(|&i| i==essential_vertex || graph.degree_in_maximal_tree(i) <= 2 )
        .filter_map(|i| points.binary_search_by_key(&i, |p| p[0] ).ok() )
        .collect::<Vec<_>>();

    // if there is no travelling points in the stem, then return
    if points_in_stem_idx.is_empty() {
        return MorsePath::new(start, end, graph);
    }

    // if there is a vacant branch, then push as many points in stem toward the branch and return.
    if let Some((&vacant_branch, _)) = graph.next_vertices[essential_vertex].iter()
        .zip( graph.next_vertices[essential_vertex].iter().skip(1).chain([next_branch].iter()) )
        .find(|&(&v, &w)| points.iter().all(|[x,_]| (v..w).contains(&x) ) ) // ToDo: change this to to binray search
    {
        let edge = [essential_vertex, vacant_branch];
        for idx in points_in_stem_idx {
            push::<true, N>(edge, idx, points);
        }
        return motions
    }

    // otherwise, we simply sort the points using 'essential_vertex'.
    let max_p = points[*points_in_stem_idx.last().unwrap()];
    let iter = graph.next_vertices[essential_vertex].iter()
        .map( |v|
            (v..).filter_map(|x|  )
        );

    let next_iters_merged = merge_sorted_dyn_by_key( next_iters, |[s,_]| *s );


    motions
}

fn collect_sorted_points<const N: usize>(essential_vertex: usize, points: &[[usize; 2]; N], graph: &GraphInformation) -> Box<dyn Iterator<Item = [usize; 2]>> {
    let points_in_stem = (0..essential_vertex).rev()
        .take_while(|v| v==&essential_vertex || graph.degree_in_maximal_tree(*v)<=2 )
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


// 'fn cmp' determines order among points during the algorithm. The order changes whether we are computing the forward path
//  or the backward path. 
fn cmp<const FORWARD: bool>(p: &[usize; 2], q: &[usize; 2]) -> Ordering {
    s
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
    const N: usize = 6;
    use crate::{ graph_collection, graphics, search::GraphInformation };

    #[test]
    fn sort_points_in_stem() -> Result<(), Box<dyn std::error::Error>> {
        let (graph, embedding, name) = graph_collection::RawGraphCollection::Grid.get();
        let graph = GraphInformation::from(&graph);

        let mut points = [[0,5], [16,20], [17,19], [20,6], [21,21], [49,49]];

        let morse_path = super::sort_points_in_stem::<N, true>(20, &mut points, &graph);

        let path = morse_path.get_geodesic();

        // draw the generated paths
        println!("drawing paths...");
        graphics::draw_edge_path::<N>(&path, &"sort_points_in_stem_test", &name, &embedding, &graph)?;

        Ok(())
    }
}

