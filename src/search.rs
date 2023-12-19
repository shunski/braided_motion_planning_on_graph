use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::ops::RangeBounds;

use topo_spaces::cubical::{EdgePath, Cube};
use alg::non_commutative::permutation::{ConstPermutation, VecPermutation};
use alg::non_commutative::Group;
use topo_spaces::graph::{CubeFactor, RawSimpleGraph};

use crate::{MorsePath, MorseCube, SortedArray};

#[derive(Debug, PartialEq, Eq, Clone)]
struct Node<const N: usize> {
    permutation: ConstPermutation<N>,
    path: EdgePath,
}

pub struct UcsHandle<const N: usize> {
    data: BinaryHeap<Node<N>>,
    atomic_paths: Vec<EdgePath>,
    considered_paths: Vec<EdgePath>,
}

impl<const N: usize> UcsHandle<N> {
    pub fn new(mut atomic_paths: Vec<EdgePath> ) -> Self {
        let mut inverses = atomic_paths.clone().into_iter().map(|p| p.inverse()).collect::<Vec<_>>();
        atomic_paths.append(&mut inverses);
        Self { data: BinaryHeap::new(), atomic_paths, considered_paths: Vec::new() }
    }
    
    pub fn search(mut self, goal: ConstPermutation<N>) -> EdgePath {
        let mut curr_node = Node{permutation: ConstPermutation::identity(), path: EdgePath::trivial_path()};
        while curr_node.permutation != goal { 
            // || curr_node.path.len() > 12 {
            // println!("search loop");
            
            // if 'curr_node' is not a path, then compute the cost for each edge (i.e. length of geodesics)
            // and add them to the heap
            for atomic_path in self.atomic_paths.iter() {
                let new_path = curr_node.path.clone().composed_with( atomic_path.clone() );
                let new_path = new_path.reduce_to_geodesic();
                if new_path.len() ==1 || self.considered_paths.contains(&new_path) {
                    continue;
                }
                let s = new_path.evaluate_permutation();

                self.data.push(Node{ permutation: s, path: new_path });
            }
            self.considered_paths.push(curr_node.path);
            curr_node = self.data.pop().unwrap();
            // println!("current path len = {:?}", curr_node.path.len());
        }
        // let list = self.data.clone().into_sorted_vec().into_iter().map( |node| node.path.len() ).collect::<Vec<_>>();
        // println!("length list = {:?}", list);

        println!("path = {:?}, \npath.len() = {}", curr_node.path, curr_node.path.len());
        curr_node.path
    }
}

impl<const N: usize> std::cmp::PartialOrd for Node<N> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(&other))
    }
}

impl<const N: usize> std::cmp::Ord for Node<N> {
    // The shorter the geodesic, the "greater" the path is in this ordering.
    // This reverse ordering is needed for min heap.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.path.len().cmp(&self.path.len())
    }
}

pub fn dynamic_search_on_tree(graph: &RawSimpleGraph, swapping: VecPermutation, start_pos: &[usize]) -> EdgePath {
    // get the path from the start to the base point.
    let start_vertex = Cube::from(start_pos.iter().map(|&v| CubeFactor::Vertex(v) ).collect::<Vec<_>>());
    let geodesic_from_start = graph.get_edge_path_to_from_base_pt(start_vertex, true);

    // get the path from the base point to the goal.
    let end_vertex = Cube::from((0..start_pos.len()).map(|i| CubeFactor::Vertex(start_pos[swapping.eval(i)]) ).collect::<Vec<_>>());
    let geodesic_to_goal = graph.get_edge_path_to_from_base_pt(end_vertex, false);

    // find the first essential vertex
    let first_essential_vertex = (0..).find(|&i| graph.degree_of(i) > 2).unwrap();

    // create the base vertex
    let base_vertex = (0..start_pos.len()).map(|i| swapping.eval(i)).collect::<Vec<_>>();
    println!("base vertex = {base_vertex:?}");

    // run the main recursive algorithm
    let cells = dynamic_search_on_tree_recc(graph, &base_vertex, first_essential_vertex, start_pos);

    println!("critical cells: {cells:?}");

    // add enough zero cells
    let cells = cells.into_iter().map(|(mut c, dir)| {
        let mut vertices_blocked_by_base = (0..start_pos.len() - c.len()).map(|i| CubeFactor::Vertex(i) ).collect::<Vec<_>>();
        vertices_blocked_by_base.append(&mut c);
        (Cube::from( vertices_blocked_by_base ), dir)
    }).collect::<Vec<_>>();

    let path = cells.into_iter().fold(geodesic_from_start, |accum, (cell, direction)| {
            let path = graph.get_edge_path(cell);
            let path = if direction {
                path
            } else {
                path.inverse()
            };

            accum.composed_with(path)
        });

    let path = path.composed_with(geodesic_to_goal);
    println!("generated path: {path:?}");
    path
}

fn dynamic_search_on_tree_recc(graph: &RawSimpleGraph, points: &[usize], essential_vertex: usize, start_pos: &[usize]) -> Vec<(Cube, bool)> {
    // create the set of vertices that comes right after the current essential vertex. 
    // the size of the set has to be equal to the number of branches - 1
    let next_vertices = graph.adjacent_vertices_iter(essential_vertex).filter(|&i| i > essential_vertex).collect::<Vec<_>>();
    println!("next_vertices = {next_vertices:?}");

    // 'first_returning_idx' is the index of the smallest point that came from another branch and goes back to the same branch
    let first_returning_idx = (0..points.len())
        .filter(|&idx| start_pos[idx] > essential_vertex && start_pos[points[idx]] > essential_vertex )
        .find(|&idx| {
            (0..next_vertices.len()-1).any(|i| 
                (next_vertices[i]..next_vertices[i+1]).contains(&start_pos[idx]) &&
                (next_vertices[i]..next_vertices[i+1]).contains(&start_pos[points[idx]])
            )
        })
        .unwrap_or(points.len());

    println!("first returning idx = {first_returning_idx}");

    // create 'swapping_orders', i.e., the swappings that we leave for other branches to do.
    let mut swapping_orders = vec![0; next_vertices.len()];
    let mut curr_next_vertices_idx = 0;
    let mut curr_branch = next_vertices[0]..next_vertices[1];
    for i in first_returning_idx..points.len() {
        while !curr_branch.contains(&start_pos[i]) {
            // if 'curr_branch' does not contain 'i', then update the 'curr_branch' to be the next branch.
            curr_next_vertices_idx += 1;
            let upper_bound = if curr_next_vertices_idx+1 < next_vertices.len() {next_vertices[curr_next_vertices_idx+1]} else { usize::MAX };
            curr_branch = next_vertices[curr_next_vertices_idx]..upper_bound;
        }
        swapping_orders[curr_next_vertices_idx] += 1;
    }
    println!("swapping orders = {swapping_orders:?}");

    // order the swapping by the dynamic programming, and swap the others on the current essential vertex
    let mut paths_from_branches = Vec::new();
    let next_essential_vertices = get_next_essential_vertices(essential_vertex, &graph);
    let mut curr_idx = points.len();
    for (i, (&n_orders, &next_essential_vertex)) in swapping_orders.iter().zip(next_essential_vertices.iter()).enumerate().rev() {
        curr_idx -= n_orders;
        let order_range = curr_idx..curr_idx+n_orders;
        println!("================ calling recursive functions ==============");
        let path = dynamic_search_on_tree_recc(graph, &points[order_range], next_essential_vertex, start_pos);
        println!("========================================================");
        println!("path = {path:?}");
        let change_of_basis = swap_chunk(essential_vertex, &swapping_orders[i..], &next_vertices[i..]);
        println!("change_of_basis = {change_of_basis:?}");
        let change_of_basis_inverse = change_of_basis.clone().into_iter().rev().map(|x| (x, false)).collect::<Vec<_>>();
        let change_of_basis = change_of_basis.into_iter().map(|x| (x, true)).collect::<Vec<_>>();
        let mut path = change_of_basis.into_iter()
            .chain( path.into_iter() )
            .chain( change_of_basis_inverse.into_iter() )
            .collect::<Vec<_>>();
        paths_from_branches.append( &mut path );
    }

    // do the necessary swapping for the returning points and traveling points.
    // these points have indces >= "first_returning_idx"
    // for n in swapping_orders {
        
    // }

    // get the ordering in the current vertex
    println!("swapping at current vertex. points = {:?}, next_vertices = {next_vertices:?}", Vec::from(&points[0..first_returning_idx]));
    let path_at_current_vertex = swap_at_current_vertex(Vec::from(points), essential_vertex, &next_vertices, &vec![0; next_vertices.len()]);

    // concatenate all the paths and evaluate the geodesics
    paths_from_branches.append(&mut path_at_current_vertex.into_iter().map(|p| (p, true)).collect::<Vec<_>>());

    paths_from_branches
}


// 'branch_idx' is the index of a vertex in 'next_vertex' to which the point is heading.
// 'next_vertices' is the labels of vertices
// 'blocked_vertices' is counting the points in the branch 
// returns some cube if it is the non-trivial group element
#[allow(unused)]
fn get_the_critical_one_cell_of(n_points: usize, essential_vertex: usize, branch_idx: usize, next_vertices: &[usize], blocked_vertices: &mut[usize]) -> Option<Cube> {
    // if the edge is not disrespecting anything then return 'None'.
    if (0..branch_idx).all(|i| blocked_vertices[i]==0) {
        blocked_vertices[branch_idx] += 1;
        return None 
    }

    // add all the vertices that are blocked by the base point.
    let n_cells_blocked_by_base = n_points-blocked_vertices.iter().map(|&x| x).sum::<usize>() - 1;
    let mut out: Vec<_> = (0..n_cells_blocked_by_base).map(|v| CubeFactor::Vertex(v)).collect();
    
    // Now add the edge
    out.push( CubeFactor::Edge( [essential_vertex, next_vertices[branch_idx]] ));

    // Finally, add all the vertices blocked by the edge
    for (&v, n) in next_vertices.iter().zip(blocked_vertices.iter_mut()) {
        if v == next_vertices[branch_idx] {
            for i in  v+1..v+1+*n { out.push(CubeFactor::Vertex(i)) }
            *n+=1;
        } else {
            for i in  v..v+*n { out.push(CubeFactor::Vertex(i)) }
        }
    }

    debug_assert!(n_points == out.len());
    Some(Cube::from(out))
}


// returns the sequence of motions (critical 1-cells) representing the swapping such that
// the 'blocked_vertices[0]' many points are disrespected
fn swap_chunk(essential_vertex: usize, blocked_vertices: &[usize], next_vertices: &[usize]) -> Vec<Cube> {
    debug_assert!(blocked_vertices.len() == next_vertices.len() );
    let mut zero_cell = (0..next_vertices.len())
        .map(|i| (0..blocked_vertices[i]).map(move |j| CubeFactor::Vertex(next_vertices[i] + j) ) )
        .flatten()
        .collect::<Vec<_>>();
    let n_disresprected_vertices = blocked_vertices[0];

    // create the motions
    let mut out = Vec::new();
    while zero_cell.len() > n_disresprected_vertices {
        let non_order_respecting_point = &mut zero_cell[n_disresprected_vertices];
        *non_order_respecting_point = CubeFactor::Edge([essential_vertex, non_order_respecting_point.vertex()]);
        out.push(Cube::from(zero_cell.clone()));
        zero_cell.remove(n_disresprected_vertices);
    }
    out
}

fn swap_at_current_vertex(mut points: Vec<usize>, essential_vertex: usize, next_vertices: &[usize], blocked_vertices: &[usize]) -> Vec<Cube> {
    // create a partition of points by indeces, so that each class in the partition is the maximal consecutive ordered sequence
    let mut prev_idx = 0;
    let mut partition = (0..points.len()).filter(|&i| i==points.len()-1 || points[i] > points[i+1] )
        .map(|i| {
            let out = prev_idx..i+1;
            prev_idx = i+1;
            out
        })
        .collect::<Vec<_>>();
    

    let mut out = Vec::new();
    // loop until 'points' is ordered
    while !(1..points.len()).all(|i| points[i-1] < points[i] ) {
        println!("swapping");
        // create the 'blocked_points'
        // 'blocked_points[i]' is a vector containing all the points blocked by the essential vertex that live in 'i'th branch
        let mut blocked_points = vec![Vec::new(); next_vertices.len()];
        let mut first_range = 0..0;
        for (range, branch_idx) in partition.iter().rev().zip( (0..next_vertices.len()).rev() ) {
            blocked_points[branch_idx] = Vec::from( &points[range.clone()] );
            
            first_range = range.clone();
        }

        println!("first_range = {first_range:?}");

        // do the non-order-respecting motions
        let mut motions = Vec::new();
        while let Some((idx, line)) = blocked_points.iter_mut()
            .enumerate()
            .filter(|(_, x)| !x.is_empty())
            .min_by_key(|(_, x)| x[0]) 
        {
            line.remove(0);

            
            let mut critical_cell = Vec::new();
            // first add the edge
            critical_cell.push(CubeFactor::Edge([essential_vertex, next_vertices[idx]]));
            // then add the blocked vertices
            for (branch_idx, n_blocked_vertices) in blocked_points.iter().map(|line| line.len() ).enumerate() {
                let mut i = 0;
                if idx == branch_idx { i += 1 };

                for _ in 0..n_blocked_vertices {
                    let v = next_vertices[branch_idx] + i;
                    critical_cell.push(CubeFactor::Vertex(v));
                    i += 1;
                }
                // add the points from the given blocked vertex
                for _ in 0..blocked_vertices[branch_idx] {
                    let v = next_vertices[branch_idx] + i;
                    critical_cell.push(CubeFactor::Vertex(v));
                    i += 1;
                }
            }

            println!("critical_cell = {critical_cell:?}");

            // check whether or not the cell is critical. if it is not critical then do not add the cell and skip
            if critical_cell.len()<2 || critical_cell[1].vertex() > next_vertices[idx] {
                continue;
            }

            
            motions.push( Cube::from(critical_cell) );
        }

        // Because some of the vertices are now ordered, we have to modify the partition
        while let Some(r) = partition.pop() {
            if r == first_range {
                break;
            }
        }
        partition.push(first_range.start..points.len());

        // Because some of the vertices are now ordered, we have to sort the portion of 'points'
        points[first_range.start..].sort();


        out.append(&mut motions);
    }

    out
}

fn get_next_essential_vertices(essential_vertex: usize, graph: &RawSimpleGraph) -> Vec<usize> {
    let mut out = Vec::new();
    for mut v in graph.adjacent_vertices_iter(essential_vertex).skip(1) {
        while graph.degree_of(v) == 2 {
            v += 1;
        }

        out.push(v);
    }
    out
}

fn get_next_vertices(essential_vertex: usize, graph: &RawSimpleGraph) -> Vec<usize> {
    graph.adjacent_vertices_iter(essential_vertex)
        .skip(1)
        .filter(|&v| graph.maximal_tree_contains( [essential_vertex, v] ))
        .collect::<Vec<_>>()
}

// fn branch_index_of(v: usize, next_vertices: &[usize]) -> usize {
//     debug_assert!(v >= next_vertices[0]);
//     let idx_op = (0..next_vertices.len()-1).find(|&i| next_vertices[i] <= v && v < next_vertices[i+1] );
//     match idx_op {
//         Some(idx) => idx,
//         None => next_vertices.len()-1,
//     }
// }

#[cfg(test)]
mod dynamic_search_test {
    use topo_spaces::cubical::Cube;
    use topo_spaces::graph::CubeFactor::{Edge, Vertex};

    use super::swap_at_current_vertex;

    #[test]
    fn swap_at_current_vertex_test() {
        // test 1
        let points = vec![2,1,0];
        let n_branches = 3;
        let essential_vertex = points.len()-1;
        let next_vertices = (0..n_branches).map(|i| (essential_vertex+1) + i * (points.len()-1) ).collect::<Vec<_>>();
        let motions = swap_at_current_vertex(points, essential_vertex, &next_vertices, &[0; 3]);

        let ans = vec![
            Cube::from(vec![Edge([2,7]), Vertex(3), Vertex(5)]),
            Cube::from(vec![Edge([2,5]), Vertex(3)])
        ];

        assert_eq!(ans, motions);


        // test 2
        let points = vec![4,6,5,2,1,3,0];
        let n_branches = 3;
        let essential_vertex = points.len()-1;
        let next_vertices = (0..n_branches).map(|i| (essential_vertex+1) + i * (points.len()-1) ).collect::<Vec<_>>();
        let motions = swap_at_current_vertex(points, essential_vertex, &next_vertices, &[0; 3]);

        // 6 is the essential vertex.
        // 7, 13, 19 are the next vertices.
        let ans = vec![
            // the function should swap the tail [2,1,3,0] first.
            Cube::from(vec![Edge([6,19]), Vertex(7), Vertex(13), Vertex(14)]),
            Cube::from(vec![Edge([6,13]), Vertex(7), Vertex(14)]),

            // next this sorts [4,6,5,0,1,2,3]
            Cube::from(vec![Edge([6,19]), Vertex(7), Vertex(8), Vertex(13), Vertex(20), Vertex(21), Vertex(22)]),
            Cube::from(vec![Edge([6,19]), Vertex(7), Vertex(8), Vertex(13), Vertex(20), Vertex(21)]),
            Cube::from(vec![Edge([6,19]), Vertex(7), Vertex(8), Vertex(13), Vertex(20)]),
            Cube::from(vec![Edge([6,19]), Vertex(7), Vertex(8), Vertex(13)]),
            Cube::from(vec![Edge([6,13]), Vertex(7)]),
        ];

        assert_eq!(ans, motions);


        // test 3 -- with unbocked vertices
        let points = vec![2,1,0];
        let n_branches = 3;
        let essential_vertex = points.len()-1;
        let next_vertices = (0..n_branches).map(|i| (essential_vertex+1) + i * (points.len()-1) ).collect::<Vec<_>>();
        let motions = swap_at_current_vertex(points, essential_vertex, &next_vertices, &[1; 3]);

        // 2 is the essential vertex.
        // 3, 5, 7 are the next vertices.
        let ans = vec![
            Cube::from(vec![Edge([2,7]), Vertex(3), Vertex(4), Vertex(5), Vertex(6), Vertex(8)]),
            Cube::from(vec![Edge([2,5]), Vertex(3), Vertex(4), Vertex(6), Vertex(7)])
        ];

        assert_eq!(ans, motions);
    }
}


pub fn path_in_tree<const N: usize>(points: [[usize; 2]; N], graph: &RawSimpleGraph) -> MorsePath<N> {
    let graph = GraphInformation::from(graph);

    let graph = ;
}

fn path_at_essential_vertex_dyn<const N: usize>(essential_vertex: usize, points: &mut [[usize; 2]; N], graph: &GraphInformation) -> MorsePath<N> {
    let mut motions = MorsePath::new();
    
    motions.compose( sort_points_in_stem::<N>(essential_vertex, points, graph) );

    for next_essential_vertex in graph.next_essential_vertices[essential_vertex] {
        motions.compose( path_at_essential_vertex_dyn(next_essential_vertex, points, graph) );
    }

    

    motions
}

fn sort_points_in_stem<const N: usize>(essential_vertex: usize, points: &mut [[usize;2]; N], graph: &GraphInformation) -> Vec<MorseCube<N>> {
    let mut motions = Vec::new();

    // 'next_branch' is the smallest vertex that is greater than any child vertices of essential vertex.
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


    let is_travelling = |idx: usize| -> bool {
        points[idx][1] <= essential_vertex || points[idx][1] >= next_branch
    };

    let cmp = |idx1: usize, idx2: usize| -> Ordering {
        if is_travelling(idx1) && !is_travelling(idx2) {
            Ordering::Less
        } else if !is_travelling(idx1) && is_travelling(idx2) {
            Ordering::Greater
        } else {
            // otherwise the ordering is determined by the goal vertex
            points[idx1][1].cmp( &points[idx2][1] )
        }
    };

    let needs_swapping = |idx: usize| -> bool {
        (0..idx).filter(|&i| !is_travelling(i))
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
        .all(|(&i, &j)| cmp(i, j) == Ordering::Less ) 
    {
        return motions;
    }


    while !points_in_stem_idx.is_empty() && !needs_swapping( *points_in_stem_idx.last().unwrap()) {
        let idx = points_in_stem_idx.pop().unwrap();
        let branch = *graph.next_vertices[essential_vertex].iter()
            .take_while(|&&i| i < points[idx][1] )
            .last().unwrap();

        // compute the motion and record it.
        motions.push( push::<true, N>([essential_vertex, branch], idx, points) );
    }

    // Having processed the trivial moves, do the swapping at the essential vertex (only once)
    let n_branches = graph.next_vertices[essential_vertex].len();
    let chunks = {
        let mut cuts = points_in_stem_idx.iter().rev().skip(1)
            .zip( points_in_stem_idx.iter().rev() )
            .filter(|(&x, &y)| cmp(x, y) == Ordering::Greater )
            .map(|(&x,_)| x )
            .take(n_branches)
            .collect::<Vec<_>>();
        cuts.reverse();
        cuts.push(usize::MAX);
        cuts.iter().zip(cuts.iter().skip(1)).map(|(&x, &y)| x..y ).collect::<Vec<_>>()
    };
    debug_assert!(chunks.len() <= n_branches);

    // push all the points to above the essential vertex
    for (&chunk, &next_vertex) in chunks.iter().rev().zip(graph.next_vertices[essential_vertex].iter().rev()) {
        for i in chunk.into_iter().rev() {
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
        cmp(idx1, idx2) 
    });

    while !falling_points.is_empty() {
        let (v, iter, next_vertex) = falling_points.first_mut().unwrap();
        let idx = points.binary_search_by_key(v, |[s,_]| *s ).unwrap();
        motions.push( push([*next_vertex, essential_vertex], idx, points));
        if let Some(w) = iter.next() {
            *v = w;

            // sort the 'falling_points' again.
            // ToDo: This sort can be more efficient.
            falling_points.sort_by(|(v,_,_), (w,_,_)| {
                let idx1 = points.binary_search_by_key(v, |[s,_]| *s ).unwrap();
                let idx2 = points.binary_search_by_key(w, |[s,_]| *s ).unwrap();
                cmp(idx1, idx2) 
            });
        } else {
            // if 'iter.next()' return 'None', then remove the tuple from 'falling_points'.
            falling_points.remove(0);
        }
    }

    motions.append( &mut sort_points_in_stem(essential_vertex, points, graph) );

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