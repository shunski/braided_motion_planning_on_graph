use braided_motion_planning::MorseCube;
use braided_motion_planning::graph_collection;
use braided_motion_planning::graphics;


// 'N' is the number of robots
const N: usize = 5;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (graph, embedding, name) = graph_collection::RawGraphCollection::Grid.get();
    
    let cube = MorseCube::<N, 1>{
        cube: [[70,105]],
        n_points_stacked_at_basepoint: 2,
        n_points_stacked_at_cube: [[vec![1,1],vec![]]],
    };

    let path = cube.get_edge_path(&graph);


    // draw the generated paths
    println!("drawing paths...");
    graphics::draw_edge_path::<N>(&path, &"test", &name, &embedding, &graph)?;

    Ok(())
}