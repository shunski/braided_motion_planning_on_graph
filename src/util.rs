use std::cmp::Ordering;

use crate::augmented_graph::AugmentedGraph;
use crate::{MorseCube, MorsePath, SortedArray};

// 'Point' represents a point (or agent, or robot) in the the maximal tree. It contains the information of 'start' and 'goal' of a point.
// This struct should only be used in the search algorithm in the maximal tree.
// This struct is in the separate module from 'search' because we want to limit the public accesses to the members.
#[derive(Clone, Copy)]
pub struct Point {
    // 'start' is the start vertex in the graph.
    start: usize,

    // 'goal' is the goal vertex in the graph.
    goal: usize,

    // 'self.meet' indicates the infimum of vertices contained in the geodesic between 'self.start' and 'self.goal' in the maximal tree.
    // 'self.meet' must be updated at the end of each constructor or each method that takes in the mutable reference.  
    meet: usize,
}



// Public Methods (Interfaces)
impl Point {
    pub fn new([start, goal]: [usize; 2], graph: &AugmentedGraph ) -> Self {
        let mut out = Self { start, goal, meet: 0 };
        out.update_meet(graph);
        out
    }

    // 'fn cmp' computes the ordering among points.
    pub fn cmp(&self, other: &Self, graph: &AugmentedGraph) -> Ordering {
        let goal_of_self = if self.meet < std::cmp::min(self.start, self.goal) {
            self.meet
        } else {
            self.goal
        };

        let goal_of_other = if other.meet < std::cmp::min(other.start, other.goal) {
            other.meet
        } else {
            other.goal
        };

        let out = goal_of_self.cmp(&goal_of_other);
        
        // 'out' can be equal if both 'goal_of_self' and 'goal_of_other' are meets. 
        // In this case we compare the starts by convention. (We could have compared the goals as well. )
        if out == Ordering::Equal {
            self.start.cmp(&other.start)
        } else {
            out
        }
    }


    // 'fn current_pos' returns the current position of the point of a graph. This depends on the direction of the algorithm,
    // i.e. 'FORWARD'
    #[inline(always)]
    pub fn curr_pos<const FORWARD: bool>(&self) -> usize {
        if FORWARD {self.start} else {self.goal}
    }
}


// This function takes a point by 'idx: usize' and its motion by 'edge: [usize; 2]' and computes the "push" motion
// in the 'UPWARD' direction, that is, 'points[idx]' will be moved from
//     (1) 'edge[0]' to 'edge[1]' if 'UPWARD' is true, and 
//     (2) 'edge[1]' to 'edge[0]' otherwise.
// If the terminal vertex (which is either 'edge[0]' or 'edge[1]') is occupied by another vertex, then it recursively
// calls this function, and try to make some space. If it is impossible to make the space, then it returns an error (i.e. Err(()) ).
// Otherwise, it will return the 'MorsePath<'_, N>'.
//
// Notes on the Inputs: This method requires that 
//     (1) 'edge' be sorted, that
//     (2a) 'FORWARD' is true imply 'points' is sorted by '_.start', and that
//     (2b) 'FORWARD' is false imply 'points' is sorted by '_.goal'.
// 'points[idx]' need not to coinside with the initial vertex (= 'edge[if UPWARD{0} else {1}]'). However,
// the geodesic between them must not contain any other point nor an essential vertex.
// The checks on these coditions will be done in the function.
// 
// Notes on the Mutable Reference 'points': When a non-order-respecting motion occurs, 'points' will be sorted by '_.start' 
// if 'FORWARD' is true and by '_.second' otherwise. No change will be applied when the output is 'Err(())'.
// 
// Notes on the Outputs: If the output of this function is 'Ok(_)', then the 'MorsePath<_,N>' containing the essential motions
// will be returned.
pub fn push<'a, const FORWARD: bool, const UPWARD: bool, const N: usize>( 
    edge: [usize; 2], 
    idx: usize, 
    points: &mut [Point; N], 
    graph: &'a AugmentedGraph
) 
    -> Result<MorsePath<'a, N>,()> 
{
    // check the input condition (1) avove.
    assert!(
        edge[0] < edge[1],
        "'edge' is not sorted."
    );

    // check the input conditions (2a) and (2b) above.
    assert!(
        points.iter().map(|p| p.curr_pos::<FORWARD>())
            .zip(points.iter().map(|p| p.curr_pos::<FORWARD>()).skip(1))
            .all(|(x,y)| x<y),
        "'points' is not sorted by {}.",
        if FORWARD {"start" } else {"goal"}    
    );


    // 'if FORWARD {points[idx].start} else {points[idx].end}' moves from 'initial' to 'terminal'.
    let [initial, terminal] = if UPWARD {edge} else {[edge[1], edge[0]]};

    // check that the geodesic between 'if FORWARD {points[idx].start} else {points[idx].end}' and 'initial' does not 
    // contain any other point nor an essential vertex. (See above)
    {
        // Make sure that the geodesic does not contain any other point.
        match points.binary_search_by_key(&initial, |p| p.curr_pos::<FORWARD>()) {
            Ok(idx2) => if idx2!=idx {panic!("geodesic contains some other points")},
            Err(idx2) => if !(idx2==idx || idx2==idx+1) {panic!("geodesic contains some other points")}
        };
        
        // Make sure that the geodesic does not contain any essential vertex.
        let p = points[idx].curr_pos::<FORWARD>();
        let result1 = graph.essential_vertices.binary_search(&p);
        let result2 = graph.essential_vertices.binary_search(&initial);
        if 
            (result1.is_ok() && result2.is_err()) || 
            (result1.is_err() && result2.is_ok()) ||
            (result1.is_ok() && result1.ok() != result2.ok() ) ||
            (result1.is_err() && result1.err() != result2.err() )
        {
            panic!("geodesic contains some essential vertex")
        }
    }


    // Create 'out', which is the output of this function
    let start = points.iter().map(|p| p.curr_pos::<FORWARD>()).collect::<Vec<_>>().try_into().unwrap();
    let out = MorsePath::new(start, start, graph);
    

    // if 'terminal' is occupied by another point, then try pushing it and make some space.
    // if it is not possible to push the point, then return 'Err(())'.
    if let Ok(i) = points.binary_search_by_key(&terminal, |p| p.curr_pos::<FORWARD>() ) {
        let mut out = if let Ok(path) = push::<FORWARD, UPWARD, N>(edge, i, points, graph) {
            path
        } else {
            return  Err(());
        };
    }


    // create the motion
    let vertices: SortedArray<N, usize> = points.iter()
        .map(|p| p.curr_pos::<FORWARD>() )
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


// private methods
impl Point {
    fn update_meet(&mut self, graph: &AugmentedGraph ) {
        // convert 'self.start' and 'self.goal' to the supremum of essential vertices less than them.
        let mut start = match graph.essential_vertices.binary_search(&self.start) {
            Ok(idx) => graph.essential_vertices[idx],
            Err(idx) => graph.essential_vertices[idx-1],
            // This does not panic because '0' is an essential vertex
        };

        let mut goal= match graph.essential_vertices.binary_search(&self.goal) {
            Ok(idx) => graph.essential_vertices[idx],
            Err(idx) => graph.essential_vertices[idx-1],
            // This does not panic because '0' is an essential vertex
        };

        // Now find the 'meet' by tracing down the parents.
        while start != goal {
            let greater_point = if start > goal {
                &mut start
            } else {
                &mut goal
            };

            *greater_point = graph.parent[*greater_point];
        }

        self.meet = self.start;
    }


    // 'fn move_by' moves 'self' by 'edge'
    // This function requires that 'edge' be sorted. This will not be checked in this method because this method is private.
    fn move_by<const UPWARD: bool>(&mut self, edge: [usize; 2], graph: &AugmentedGraph ) {
        debug_assert!(
            edge[0] < edge[1],
            "'edge' is not sorted. Check the public methods that use this method."
        );

        self.update_meet(graph);
    }
}
