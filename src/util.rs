use std::{
    cmp,
    fmt,
    mem
};

use crate::augmented_graph::AugmentedGraph;
use crate::{MorseCube, MorsePath};




// 'SortedArray' is the sorted array type that has a fixed max size ('N').
// This type does not use heap allocation, so it is copiable.
#[derive(Copy)]
pub struct SortedArray<const N: usize, T: Ord> {
    data: [T; N],
    len: usize,
}

impl<const N: usize, T: Ord + Copy + fmt::Debug> SortedArray<N, T> {}

impl<const N: usize, T: Ord + Copy> std::ops::Deref for SortedArray<N, T> {
    type Target = [T];
    fn deref(&self) -> &Self::Target {
        &self.data[..self.len]
    }
}

impl<const N: usize, T: Ord + Copy> Clone for SortedArray<N,T> {
    fn clone(&self) -> Self {
        let mut out = Self::new();
        for &x in self.iter() {
            out.push(x);
        }
        out
    }
}

impl<const N: usize, T: Ord + Copy + fmt::Debug> fmt::Debug for SortedArray<N, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.data[..self.len].fmt(f)
    }
}

impl<const N: usize, T: Ord> cmp::PartialEq for SortedArray<N, T> {
    fn eq(&self, other: &Self) -> bool {
        (&*self).eq(&*other)
    }
}

impl<const N: usize, const M: usize, T: Ord + Copy + fmt::Debug> From<[T; M]> for SortedArray<N, T> {
    fn from(mut value: [T; M]) -> Self {
        assert!( M <= N, "cannot build an array of size {N} from an array of size {M}");

        assert!(
            value.iter()
                .zip(value.iter().skip(1))
                .all(|(a,b)| a <= b ),
            "the input array is not sorted."
        );

        let mut out = Self::new();
        for x in value {
            out.push(x);
        }
        out
    }
}


impl<const N: usize, T: Ord + Copy + fmt::Debug> TryFrom<Vec<T>> for SortedArray<N, T> {
    type Error = Self;
    // if 'value.len() > N', then this function creates an array of 'N' elements that contains first 'N' elements of 'value',
    // and encloses it in 'Err'
    // This method will also panic is the iterator is not sorted.
    fn try_from(mut value: Vec<T>) -> Result<Self, Self::Error> {
        value.sort();
        let is_error = value.len() > N;
        let mut out = Self::new();
        for x in value.into_iter().take(N) {
            out.push(x);
        }
        if is_error {
            Err(out)
        } else {
            Ok(out)
        }
    }
}

impl<const N: usize, T: Ord + Copy> SortedArray<N,T> {
    #[inline(always)]
    pub fn new() -> Self {
        let data = unsafe {
            mem::MaybeUninit::uninit().assume_init()
        };
        Self { data, len: 0 }
    }


    // 'push' adds input argument at the end. If the input argument is not greater than the
    // largest element of the current cell, then this function panics. If the user is unsure whether or not this
    // condition is satisfied, then 'insert' should be used instead.
    #[inline]
    pub fn push(&mut self, v: T) {
        // check the size of the array.
        assert!(self.len<N, "The array is full.");

        // check that the input argument is greater than the largest element.
        assert!(
            self.is_empty() || self.last().unwrap() <= &v,
            "The value cannot be added because it is less than the largest element."
        );

        // The following unsafe block is safe, because we have checked that the size of the array is less than 'N'.
        unsafe{
            self.data.as_mut_ptr().add(self.len).write(v);
        }

        self.len+=1;


        debug_assert!(
            self.iter().zip(self.iter().skip(1)).all(|(x,y)| x<=y ),
            "self is not sorted."
        );
    }


    // 'fn insert' inserts the input argument into the array so that the array is still sorted.
    pub fn insert(&mut self, v: T) {
        // check the size of the array.
        assert!(self.len<N, "The array is full.");

        // The following unsafe block is safe, because we have checked that the size of the array is less than 'N'.
        unsafe{
            self.data.as_mut_ptr().add(self.len).write(v);
        }

        self.len+=1;
        
        // Sort the array.
        // ToDo: This can be more efficient.
        self.data.sort();
    }


    // 'fn remove' is analogous to that of 'Vec', i.e., it removes the element at the specified index, and returns the value.
    pub fn remove(&mut self, idx: usize) -> T {
        assert!(idx < self.len, "index out of bounds: array is of length={}, but 'idx'={idx}", self.len);

        let out = self[idx];
        for i in idx..self.len-1 {
            self[i] = self[i+1]; 
        }

        self.len -= 1;
        out
    }


    // 'fn into_iter' is the ad-hoc implementaion of into_iter.
    // ToDo: implement the IntoIter for this type (without using 'into_iter' of 'data').
    pub fn into_iter(self) -> impl Iterator<Item = T> { 
        self.data.into_iter().take(self.len)
    }
}



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
    pub fn cmp(&self, other: &Self, graph: &AugmentedGraph) -> cmp::Ordering {
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
        if out == cmp::Ordering::Equal {
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
// The checks on these conditions will be done in the function.
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


    // Create 'out', which is the output of this function.
    // This will record the movements done in this function.
    let start = points.iter().map(|p| p.curr_pos::<FORWARD>()).collect::<Vec<_>>().try_into().unwrap();
    let mut out = MorsePath::new(start, start, graph);
    

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
        .filter(|&s| s!=points[idx].curr_pos())
        .collect::<Vec<_>>()
        .try_into().unwrap();
    debug_assert!(vertices.len() == N-1);

    let edges: SortedArray<N, [usize; 2]> = [edge].into();


    // Create the motion and record it if it is critical.
    let motion = MorseCube::new_unchecked(edges, vertices);
    if motion.flows_to_critical_one_cell(&graph) {
        out.add( motion );
    }


    // finally, physically move the point
    points[idx].move_by::<UPWARD>(edge, graph);
    points.sort_by_key(|p| p.curr_pos::<FORWARD>()); // ToDo: This can be more efficient

    Ok(out)
}


// private methods
impl Point {
    // 'fn update_meet' updates the "meet" of 'start' and 'goal'.
    // This private method is typically used by other method when the struct 'Point' 
    // is constructed or its 'start' and 'goal' values are modified.
    // This method changes the value of 'meet' only, and the values of 'start' and 'goal' remain unchanged.
    fn update_meet(&mut self, graph: &AugmentedGraph ) {
        // update the meet.
        self.meet = graph.meet_of(self.start, self.goal);
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
