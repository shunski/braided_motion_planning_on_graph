use std::collections::HashMap;

use alg::lin_alg::ConstVector;
use topo_spaces::graph;
use topo_spaces::graph::Graph;

use std::f64::consts::PI;

pub enum Collection {
    Star,
    FourTriangles,
    ATree,
    Grid,
    ThreeVerticesOfValencyThree
}

pub fn get(c: Collection) -> (String, HashMap<String, ConstVector<f64, 2>>, Graph<i64>) {
    match c {
        Collection::Star => {
            // A star graph with 'n' vertices, the star graph with three branches
            let graph = graph! {
                ["o", "b0"],
                ["o", "b1"],
                ["o", "b2"]
            };

            // graph embedding
            let graph_embedding = {
                let mut e = HashMap::new();
                e.insert( 'o'.to_string(), ConstVector::from([ 0.0, 0.0]));
                for i in 0..3 {
                    let theta = 2.0 * PI / 3.0 * i as f64;
                    e.insert( format!("b{}", i), ConstVector::from([ -theta.sin(), theta.cos() ]));
                }
                e
            };
            ("Star(3)".to_string(), graph_embedding, graph)
        },

        Collection::FourTriangles => {
            let graph = graph! {
                ["A", "c"],
                ["A", "b"],
                ["B", "c"],
                ["B", "a"],
                ["C", "a"],
                ["C", "b"],
                ["a", "c"],
                ["a", "b"],
                ["b", "c"]
            };
        
            // graph embedding
            let v1 = ConstVector::from([1.0, 0.0]);
            let v2 = ConstVector::from([0.5, 3_f64.sqrt()/2.0]);
            let o = ConstVector::from([-1.0, -0.9]);
        
        
            let graph_embedding = {
                let mut e = HashMap::new();
                e.insert( 'A'.to_string(), o+v2+v2);
                e.insert( 'B'.to_string(), o);
                e.insert( 'C'.to_string(), o+v1+v1);
                e.insert( 'a'.to_string(), o+v1);
                e.insert( 'b'.to_string(), o+v1+v2);
                e.insert( 'c'.to_string(), o+v2);
                e
            };
            ("FourTriangles".to_string(), graph_embedding, graph)
        },

        Collection::ATree => {
            // A connected tree with four vertices of valency 3
            let graph = graph! {
                ["o", "x"],
                ["o", "y"],
                ["o", "z"],
                ["x", "x1"],
                ["x", "x2"],
                ["y", "y1"],
                ["y", "y2"],
                ["z", "z1"],
                ["z", "z2"]
            };

            // graph embedding
            let v = |dir: usize| -> ConstVector<f64, 2> { 
                ConstVector::from([ (dir as f64/3.0*PI).sin(), (dir as f64/3.0*PI).cos(), ]) * 0.6
            };



            let graph_embedding = {
                let mut e = HashMap::new();
                e.insert( 'o'.to_string(), ConstVector::from([ 0.0, 0.0]));
                e.insert( 'x'.to_string(), v(0));
                e.insert( 'y'.to_string(), v(2));
                e.insert( 'z'.to_string(), v(4));

                e.insert( "x1".to_string(), v(0)+v(1) );
                e.insert( "x2".to_string(), v(0)+v(5));

                e.insert( "y1".to_string(), v(2)+v(1));
                e.insert( "y2".to_string(), v(2)+v(3));

                e.insert( "z1".to_string(), v(4)+v(3));
                e.insert( "z2".to_string(), v(4)+v(5));
                e
            };
            ("ATree".to_string(), graph_embedding, graph)
        },

        Collection::Grid => {
            // A small grid-like graph
            let graph = graph! {
                ["1", "2"],
                ["2", "3"],
                ["3", "4"],
                ["1", "4"],
                ["1", "5"],
                ["1", "6"],
                ["2", "7"],
                ["2", "8"],
                ["3", "9"],
                ["3", "10"],
                ["4", "11"],
                ["4", "12"],
                ["5", "6"],
                ["6", "7"],
                ["7", "8"],
                ["8", "9"],
                ["9", "10"],
                ["10", "11"],
                ["11", "12"],
                ["5", "12"]
            };

            // graph embedding
            let graph_embedding = {
                let mut e = HashMap::new();
                let e1 = ConstVector::from([ 0.6, 0.0]);
                let e2 = ConstVector::from([ 0.0, 0.6]);
                let o = ConstVector::from([ -0.9, -0.9]);
                e.insert(  "1".to_string(), o + e1*2.0 + e2*2.0 );
                e.insert(  "2".to_string(), o + e1*1.0 + e2*2.0 );
                e.insert(  "3".to_string(), o + e1*1.0 + e2*1.0 );
                e.insert(  "4".to_string(), o + e1*2.0 + e2*1.0 );
                e.insert(  "5".to_string(), o + e1*3.0 + e2*2.0 );
                e.insert(  "6".to_string(), o + e1*2.0 + e2*3.0 );
                e.insert(  "7".to_string(), o + e1*1.0 + e2*3.0 );
                e.insert(  "8".to_string(), o + e1*0.0 + e2*2.0 );
                e.insert(  "9".to_string(), o + e1*0.0 + e2*1.0 );
                e.insert( "10".to_string(), o + e1*1.0 + e2*0.0 );
                e.insert( "11".to_string(), o + e1*2.0 + e2*0.0 );
                e.insert( "12".to_string(), o + e1*3.0 + e2*1.0 );
                e
            };
            ("Grid".to_string(), graph_embedding, graph)
        }

        Collection::ThreeVerticesOfValencyThree => {
            // A connected tree with three vertices of valency 3
            let graph = graph! {
                ["o1", "o2"],
                ["o2", "b1"],
                ["o2", "b2"],
                ["a1", "b1"],
                ["b1", "c1"],
                ["a2", "b2"],
                ["b2", "c2"]
            };

            // graph embedding
            let graph_embedding = {
                let mut e = HashMap::new();
                let e1 = ConstVector::from([ 0.9, 0.0]);
                let e2 = ConstVector::from([ 0.0, 0.9]);
                let o = ConstVector::from([ 0.0, -0.9]);
                e.insert(  "o1".to_string(), o );
                e.insert(  "o2".to_string(), o + e2 );
                e.insert(  "a1".to_string(), o + e1 );
                e.insert(  "b1".to_string(), o + e1 + e2 );
                e.insert(  "c1".to_string(), o + e1 + e2 * 2.0 );
                e.insert(  "a2".to_string(), o - e1 );
                e.insert(  "b2".to_string(), o - e1 + e2 );
                e.insert(  "c2".to_string(), o - e1 + e2 * 2.0 );
                e
            };
            ("ThreeVerticesOfValencyThree".to_string(), graph_embedding, graph)
        }
    }
}