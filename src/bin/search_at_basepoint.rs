use std::time::Instant;

use alg::non_commutative::permutation::ConstPermutation;
use topo_spaces::cubical::UdnMorseCplx;

use braided_motion_planning::{
    graph_collection,
    graphics,
    search,
};

const N_ROBOTS: usize = 7;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    // constructing the base graph and setting the goal
    let (name, embedding, graph) = graph_collection::get(graph_collection::Collection::ATree);
    let goal = ConstPermutation::from([1,4,0,5,3,2,6]);
    
    // starting the timer for search algorithm
    let search_start = Instant::now();

    // constructing the cplx
    let cplx = UdnMorseCplx::new(&graph, N_ROBOTS, Some(1));

    // modify the embedding
    let modified_embedding = cplx.factor_pool.modify_embedding(embedding);
    
    // get the paths
    let atomic_paths = cplx.get_cells_of_dim(1)
        .into_iter()
        .map(|c| cplx.factor_pool.get_edge_path(c) )
        .collect::<Vec<_>>();

    // run the search
    let handle = search::UcsHandle::new(atomic_paths);
    let geodesic = handle.search(goal);


    // print the time elapsed to complete the search
    let search_time = search_start.elapsed();
    println!("The search took {:.3} seconds. drawing paths...", search_time.as_secs_f64());

    // draw the generated paths
    graphics::draw_paths::<N_ROBOTS>(&[geodesic], &modified_embedding, &cplx, &name, "SearchAtBasePoint")?;

    Ok(())
}
