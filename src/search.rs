use std::cmp::Ordering;

use quasi_iter::merge_sorted::{merge_sorted_dyn_by_key, merge_sorted_dyn_by, MergeSortedDynBy};
use topo_spaces::graph::RawSimpleGraph;
use crate::augmented_graph::AugmentedGraph;
use crate::{MorsePath, MorseCube, SortedArray};


// The main algorithm that computes the path in the maximal tree.
pub fn path_in_tree<'a, const N: usize>(points: &mut [[usize; 2]; N], graph: &'a RawSimpleGraph) -> MorsePath<'a, N> {
    let graph = AugmentedGraph::<'a>::from(graph);

    // the dynamic programming and recursively sort the points in the graph at each essential vertex.
    let first_vertex_of_deg_greater_than_two = *graph.next_essential_vertices[0].first().unwrap();
    let mut motions = path_at_essential_vertex_dyn(first_vertex_of_deg_greater_than_two, points, &graph);

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
    essential_vertex: usize, points: &mut [[usize; 2]; N], graph: &'b AugmentedGraph<'a>
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

fn sort_points_in_stem<'a, const N: usize, const FORWARD: bool>(essential_vertex: usize, points: &mut [[usize;2]; N], graph: &'a AugmentedGraph) -> MorsePath<'a, N> {
    let (start, end): (Vec<_>, Vec<_>) = points.iter().copied().map(|[x,y]| (x,y) ).unzip();
    let start = start.try_into().unwrap();
    let end = end.try_into().unwrap();
    let mut motions: MorsePath<'_, N> = MorsePath::new(start, end, graph);

    // 'next_branch' is the smallest vertex that is greater than any child vertices of the essential vertex.
    let next_branch = *graph.next_vertices[graph.parent[essential_vertex]].iter()
        .find(|&&v| v>essential_vertex)
        .unwrap_or(&usize::MAX);


    // this closure returns 'true' if the swapping has to be done at the current essential vertex.
    // That is, if this closure returns 'false', then we can simply push this vertex to the next child branch.
    let needs_swapping = |idx: usize, points: &[[usize;2];N]| -> bool {
        (0..idx).filter(|&i| local_cmp(i, idx, points, essential_vertex, graph).is_gt())
            .any(|i|
                // the following expression returns ture if 
                // (1) 'points[idx]' is travelling (this automatically implies 'points[i]' is also travelling), or if
                is_travelling(idx, points, essential_vertex, graph)
                ||
                // (2) the following expression returns true if the goal of 'i' and the goal 'idx' belong to the same child branch
                graph.next_vertices[essential_vertex]
                    .iter()
                    .zip( graph.next_vertices[essential_vertex].iter().skip(1).chain(&[next_branch])) 
                    .any( |(&a, &b)| (a..b).contains(&points[i][1]) && (a..b).contains(&points[idx][1]) )
            )
    };


    // 'points_in_stem_idx' is the points located at vertices less than or equal to the current essential vertex and 
    // larger than the previous essential vertex of valency greater than 3 or greater than or equal to the previous
    // essential vertex of valency 1.
    let mut points_in_stem_idx = (0..=essential_vertex).rev()
        .take_while(|&i| i==essential_vertex || graph.degree_in_maximal_tree(i) <= 2 )
        .filter_map(|i| points.binary_search_by_key(&i, |p| p[0] ).ok() )
        .collect::<Vec<_>>();
    points_in_stem_idx.reverse();
    

    println!("points_in_stem_idx={points_in_stem_idx:?}");

    while !points_in_stem_idx.is_empty() && !needs_swapping( *points_in_stem_idx.last().unwrap(), points) && !is_travelling( *points_in_stem_idx.last().unwrap(), points, essential_vertex, graph) {
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
        .all(|(&i, &j)| local_cmp(i, j, points, essential_vertex, graph).is_lt() ) 
        &&
        points_in_stem_idx.iter().all(|&i| is_travelling(i, points, essential_vertex, graph) )
    {
        return motions;
    }

    // Having processed the trivial moves, do the swapping at the essential vertex (only once)
    let n_branches = graph.next_vertices[essential_vertex].len();
    let swapping_chunks = {
        let mut cuts = points_in_stem_idx.iter().rev().skip(1)
            .zip( points_in_stem_idx.iter().rev() )
            .filter(|(&x, &y)| local_cmp(x, y, points, essential_vertex, graph).is_gt() )
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
        swapping_chunks.clone().into_iter().all(|chunk| chunk.end == points_in_stem_idx.len() || local_cmp(chunk.end-1, chunk.end, points, essential_vertex, graph).is_gt() ),
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
                local_cmp(idx1, idx2, &curr_points, essential_vertex, graph)
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


// 'fn sort_travelling_points' sorts points in the stem, and send those points that need to be sorted
// in the child branches to the child branches.
// This function requires that 
//     (1) each point in 'points' inside the stem needs to be a travelling point,
//     (2) and that 'points' is sorted by the first (or second) element if 'FORWARD'='true' (or ='false').
fn sort_travelling_points<'a, const N: usize, const FORWARD: bool>( essential_vertex: usize, points: &mut [[usize; 2]; N], graph: &'a AugmentedGraph ) -> MorsePath<'a, N> {
    // First, check the two conditions specified above.
    assert!(
        (0..points.len()).all(|i| is_travelling::<N,FORWARD>(i, points, essential_vertex, graph)),
        "Not all the points are travelling."
    );
    assert!(
        points.iter()
            .map(|p| if FORWARD {p[0]} else {p[1]} )
            .zip(points.iter().skip(1).map(|p| if FORWARD {p[0]} else {p[1]} ))
            .all(|(x,y)| x < y),
        "Points not sorted with repspect to the {} element.",
        if FORWARD {"first"} else {"second"}
    );


    // initialize 'motions', which is the output of this function.
    let start = points.iter().map(|&[s,_]| s).collect::<Vec<_>>().try_into().unwrap();
    let end = points.iter().map(|&[_,g]| g).collect::<Vec<_>>().try_into().unwrap();
    let mut motions = MorsePath::new(start, end, &graph);
    
    // 'next_branch' is the smallest vertex that is greater than any child vertices of the essential vertex.
    let next_branch = *graph.next_vertices[graph.parent[essential_vertex]].iter()
        .find(|&&v| v>essential_vertex)
        .unwrap_or(&usize::MAX);

    
    // collect the points in stem.
    let points_in_stem_idx = (0..=essential_vertex).rev()
        .take_while(|&i| i==essential_vertex || graph.degree_in_maximal_tree(i) <= 2 )
        .filter_map(|i| points.binary_search_by_key(&i, |p| p[0] ).ok() )
        .collect::<Vec<_>>();

    // if there is no travelling points in the stem, then return
    if points_in_stem_idx.is_empty() {
        return MorsePath::new(start, end, graph);
    }

    // if there is a vacant branch, then push as many points in stem toward that branch and return.
    if let Some((&vacant_branch, _)) = graph.next_vertices[essential_vertex].iter()
        .zip( graph.next_vertices[essential_vertex].iter().skip(1).chain([next_branch].iter()) )
        .find(|&(&v, &w)| points.iter().all(|[x,_]| !(v..w).contains(&x) ) ) // ToDo: change this to to binray search
    {
        let edge = [essential_vertex, vacant_branch];
        for idx in points_in_stem_idx {
            motions.add( push::<true, N>(edge, idx, points) );
        }
        return motions
    }

    // if there is a branch such that every point in it is greater than the largest travelling point (in the local order), 
    // then push all the points towards the branch and return.
    if let Some((&nice_branch, _)) = graph.next_vertices[essential_vertex].iter()
        .zip( graph.next_vertices[essential_vertex].iter().skip(1).chain([next_branch].iter()) )
        .find(|&(&v, &w)| points.iter().all(|[x,_]| 
            !(v..w).contains(&x) ||
            local_cmp(*x, *points_in_stem_idx.last().unwrap(), points, essential_vertex, graph).is_gt() ) 
        )
    {
        let edge = [essential_vertex, nice_branch];
        for idx in points_in_stem_idx {
            motions.add( push::<true, N>(edge, idx, points) );
        }
        return motions
    }


    // Otherwise, we have to push all the travelling points to the child branches, and leave them to swap the points.
    // The straightforward implementation is to equally distribute points to the branches with more than one essential vertices.
    for (&idx, &terminal) in points_in_stem_idx.iter().rev()
        .zip(
            graph.next_vertices[essential_vertex].iter()
                // make sure that we filter out the branches homeomorphic to a line segment
                .filter(|&&v| (v..).find(|&w| graph.degree_in_maximal_tree(w)!=2).unwrap() > 2 )
        )
        .cycle()
    {
        let initial = points[idx][if FORWARD {0} else {1}];
        motions.add(push::<true, N>([initial, terminal], idx, points))
    }


    motions
}

// This function returns an iterator over all the elements of 'points' in the order determined by their goals.
// This function requires that
//     (1) each element of 'points' be either travelling or staying in the same stem, and that
//     (2) elements of 'points' in each branch is sorted by 'local_cmp'.
fn collect_sorted_points<const N: usize>(essential_vertex: usize, points: &[[usize; 2]; N], graph: &AugmentedGraph) -> Box<dyn Iterator<Item = [usize; 2]>> {
    // check condition (1)
    assert!(
        points.iter().map(|p|  ),
        ""
    );

    let points_in_stem = (0..=essential_vertex).rev()
        .take_while(|&v| v==essential_vertex || graph.degree_in_maximal_tree(v)<=2 )
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


// 'fn is_travelling' determines whether or not the point 'p' needs to travel towards other branches through the given 'essential_vertex'.
// Note that if a point stays in the stem, then this point is by definition NOT a travelling point.
fn is_travelling<const N: usize, const FORWARD: bool>(idx: usize, points: &[[usize;2]; N], essential_vertex: usize, graph: &AugmentedGraph) -> bool {
    // 'next_branch' is the smallest vertex greater than 'essential_vertex' that is adjacent to the parent of 'essential_vertex'.
    let next_branch = *graph.next_vertices[graph.parent[essential_vertex]].iter().find(|&&v| v>essential_vertex).unwrap_or(&usize::MAX);
    let smallest_v_in_stem = (0..essential_vertex).rev().filter(|&v| graph.degree_in_maximal_tree(v)==2).last().unwrap();
    let p = points[idx][if FORWARD {1} else {0}];
    p >= next_branch || p < smallest_v_in_stem
}


#[cfg(test)]
mod path_generation_on_tree_tests {
    // 'N' is the number of robots
    const N: usize = 6;
    use crate::{ graph_collection, graphics, search::AugmentedGraph };

    #[test]
    fn sort_points_in_stem() -> Result<(), Box<dyn std::error::Error>> {
        let (graph, embedding, name) = graph_collection::RawGraphCollection::Grid.get();
        let graph = AugmentedGraph::from(&graph);

        let mut points = [[0,5], [16,20], [17,19], [20,6], [21,21], [49,49]];

        let morse_path = super::sort_points_in_stem::<N, true>(20, &mut points, &graph);

        let path = morse_path.get_geodesic();

        // draw the generated paths
        println!("drawing paths...");
        graphics::draw_edge_path::<N>(&path, &"sort_points_in_stem_test", &name, &embedding, &graph)?;

        Ok(())
    }
}

