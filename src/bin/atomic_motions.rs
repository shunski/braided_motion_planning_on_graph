use topo_spaces::cubical::UdnMorseCplx;

use braided_motion_planning::graph_collection;
use braided_motion_planning::graphics;

const N_ROBOTS: usize = 4;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (name, embedding, graph) = graph_collection::get(graph_collection::Collection::Grid);

    let cplx = UdnMorseCplx::new(&graph, N_ROBOTS, Some(1));

    // modify the embedding
    let modified_embedding = cplx.factor_pool.modify_embedding(embedding);
    
    // // -------------------- debug ----------------------
    // let cells = cplx.get_cells_of_dim(1);
    // for cell in cells.into_iter(){
    //     println!("{cell}");
    // }
    // // -------------------------------------------------
    
    // get the paths
    let edge_paths = cplx.get_cells_of_dim(1)
        .into_iter()
        .map(|c| cplx.factor_pool.get_edge_path(c) )
        .collect::<Vec<_>>();


    // draw the generated paths
    println!("drawing paths...");
    graphics::draw_paths::<N_ROBOTS>(&edge_paths, &modified_embedding, &cplx, &name, "AtomicMotions")?;

    Ok(())
}
