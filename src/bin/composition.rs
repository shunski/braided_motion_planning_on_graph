use braided_motion_planning::MorseCube;
use braided_motion_planning::graph_collection;
use braided_motion_planning::graphics;


// 'N' is the number of robots
const N: usize = 5;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (graph, embedding, name) = graph_collection::RawGraphCollection::Grid.get();
    
    let cube1 = MorseCube::<N, 1>{
        cube: [[70,105]],
        n_points_stacked_at_basepoint: 2,
        n_points_stacked_at_cube: [[vec![1,1],vec![]]],
    };

    let cube2 = MorseCube::<N, 1>{
        cube: [[70,105]],
        n_points_stacked_at_basepoint: 2,
        n_points_stacked_at_cube: [[vec![1,1],vec![]]],
    };

    let path1 = cube1.get_edge_path(&graph);
    let path2 = cube2.get_edge_path(&graph);
    println!("path1 = {path1:?}");

    let composition = path1.composed_with(path2);

    println!("compostion = {composition:?}");

    // draw the generated paths
    println!("drawing paths...");
    graphics::draw_edge_path::<N>(&composition, &"test", &name, &embedding, &graph)?;

    Ok(())
}