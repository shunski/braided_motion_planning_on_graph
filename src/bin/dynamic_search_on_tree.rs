use std::time::Instant;
use alg::non_commutative::permutation::VecPermutation;
use braided_motion_planning::{
    search::dynamic_search_on_tree,
    graph_collection,
    graphics
};

use topo_spaces::cubical::UdnMorseCplx;

const N_ROBOTS: usize = 3;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (name, embedding, graph) = graph_collection::get(graph_collection::Collection::ThreeVerticesOfValencyThree);

    
    // modify the embedding
    let (graph, modified_embedding) = {
        let cplx = UdnMorseCplx::new(&graph, N_ROBOTS, Some(1));
        let modified_embedding = cplx.factor_pool.modify_embedding(embedding);
        (cplx.factor_pool, modified_embedding)
    };
    
    // set up the initial condition
    let pos = [6, 8, 12];
    let swapping = VecPermutation::from([1,2,0]);
    
    // starting the timer for search algorithm
    let search_start = Instant::now();
   
    // compute the motion plan
    let path = dynamic_search_on_tree(&graph, swapping.clone(), &pos);
    let path = path.reduce_to_geodesic();
   
    // print the time elapsed to complete the search
    let search_time = search_start.elapsed();
    println!("The search took {:.3} seconds. drawing paths...", search_time.as_secs_f64());

    // draw the generated paths
    println!("drawing the path...");
    let start = pos;
    let end = [pos[swapping.eval(0)], pos[swapping.eval(1)], pos[swapping.eval(2)]];
    graphics::draw_path_with_endpoints::<N_ROBOTS>(&path, &start, &end, &modified_embedding, &graph, &name, "DynamicProgrammingOnTrees")?;
    Ok(())
}