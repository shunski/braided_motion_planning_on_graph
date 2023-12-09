use topo_spaces::cubical::{self, Cube, UdnMorseCplx};
use topo_spaces::graph::{CubeFactor, SimpleVertex, SimpleGraph, RawSimpleGraph};
use alg::lin_alg::ConstVector;
use crate::CubicPath;

use std::collections::HashMap;

use plotters::prelude::*;

const EDGE_PATH_ANIMATION_RESOLUTION: usize = 10;
const FLAME_DELAY: u32 = 5;

struct PathRenderer<'a, const N: usize> {
    iter: std::slice::Iter<'a, Vec<[usize; 2]>>,

    // current edge cube
    // [[start_x, start_y], [end_x, end_y]]
    curr_edge: [[ConstVector<f64,N>; 2]; 2],

    // terminal vertex of the current edge
    terminal_of_curr_edge: [usize; N],

    // current position in the current edge cube. Must satisfy 0 <= 'curr_pos' < 'EDGE_PATH_ANIMATION_RESOLUTION' at any time.
    curr_pos: usize,

    // graph embedding: a map from the vertices of the graph to the plane
    graph_embedding: &'a HashMap<SimpleVertex, ConstVector<f64, 2>>,

    vertices_labeling: &'a [SimpleVertex],
}

impl<'a, const N: usize> PathRenderer<'a, N>
{
    fn new(path: &'a Vec<Vec<[usize; 2]>>, start: Option<&[usize]>, graph_embedding: &'a HashMap<SimpleVertex, ConstVector<f64, 2>>, vertices_labeling: &'a [SimpleVertex] ) -> Self {
        let mut iter = path.iter();
        let initial_path = iter.next();

        let start = if let Some(start) = start {
            Cube::from( start.iter().map(|&v| CubeFactor::Vertex(v)).collect::<Vec<_>>() )
        } else {
            // if the 'start' is not specified, start from the base point.
            Cube::from((0..N).map(|v| CubeFactor::Vertex(v)).collect::<Vec<_>>())
        };

        let terminal_of_curr_edge = {
            let mut terminal_of_curr_edge = [0; N];
            let edges = initial_path.unwrap();
            for (i, v) in start.iter().map(|f| f.vertex()).enumerate() {
                terminal_of_curr_edge[i] = if let Some(e) = edges.iter().find(|e| e[0] == v) {
                    e[1]
                } else {
                    v
                }
            }
            terminal_of_curr_edge
        };

        let mut curr_edge = [[ConstVector::<f64, N>::zero();2];2];
        for i in 0..N {
            curr_edge[0][0][i] = graph_embedding[&vertices_labeling[start[i].vertex()]][0];
            curr_edge[0][1][i] = graph_embedding[&vertices_labeling[start[i].vertex()]][1];

            curr_edge[1][0][i] = graph_embedding[&vertices_labeling[terminal_of_curr_edge[i]]][0];
            curr_edge[1][1][i] = graph_embedding[&vertices_labeling[terminal_of_curr_edge[i]]][1];
        }

        let curr_pos = 0;
        Self {
            iter,
            curr_edge,
            terminal_of_curr_edge,
            curr_pos,
            graph_embedding,
            vertices_labeling
        }
    }
}

impl<'a, const N: usize> Iterator for PathRenderer<'a, N> {
    type Item = [ConstVector<f64, 2>; N];
    fn next(&mut self) -> Option<Self::Item> {
        // if this is the end of an edge...
        if self.curr_pos == EDGE_PATH_ANIMATION_RESOLUTION {
            // and if there is still another edge to traverse, ...
            if let Some(next_edges) = self.iter.next() {
                
                self.terminal_of_curr_edge = {
                    let mut terminal_of_next_edge = [0; N];
                    for (i, &v) in self.terminal_of_curr_edge.iter().enumerate() {
                        terminal_of_next_edge[i] = if let Some(e) = next_edges.iter().find(|e| e[0] == v) {
                            e[1]
                        } else {
                            v
                        }
                    }
                    terminal_of_next_edge
                };
                

                self.curr_edge.swap(0, 1);
                // now 'self.curr_edge[1]' has the start location of old edge.
                // rewrite the end of the new edge here.
                for i in 0..N {
                    self.curr_edge[1][0][i] = self.graph_embedding[&self.vertices_labeling[self.terminal_of_curr_edge[i]]][0];
                    self.curr_edge[1][1][i] = self.graph_embedding[&self.vertices_labeling[self.terminal_of_curr_edge[i]]][1];
                }

                self.curr_pos = 0;
            } else {
                return None;
            }
        } else {
            self.curr_pos += 1;
        }

        
        let x = (self.curr_edge[1][0] - self.curr_edge[0][0]) * self.curr_pos as f64 / EDGE_PATH_ANIMATION_RESOLUTION as f64 + self.curr_edge[0][0];
        let y= (self.curr_edge[1][1] - self.curr_edge[0][1]) * self.curr_pos as f64 / EDGE_PATH_ANIMATION_RESOLUTION as f64 + self.curr_edge[0][1];
        
        let mut out: [ConstVector<f64, 2>; N] = [ConstVector::zero(); N];
        for i in 0..N {
            out[i][0] = x[i];
            out[i][1] = y[i];
        }

        Some(out)
    }
}


pub fn draw_paths<const N: usize>(edge_paths: &[cubical::EdgePath], embedding: &HashMap<SimpleVertex, ConstVector<f64, 2>>, cplx: &UdnMorseCplx, graph_name: &str, bin_name: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut path_renderers: Vec<_> = edge_paths.iter()
        .map(|p| PathRenderer::<N>::new(
            p, 
            None, // start the motion from the base point
            &embedding, 
            cplx.factor_pool.vertices_list()
        ))
        .collect();

    let n_cols = (path_renderers.len() as f64 - 0.1).sqrt() as usize + 1;
    let n_rows = (path_renderers.len()-1) / n_cols + 1;

    let vis_root = BitMapBackend::gif(
        &format!("program_outputs/{bin_name}==={graph_name}==={N}_points.gif"), 
        (150 * n_cols as u32, 150 * n_rows as u32), 
        10
    )?.into_drawing_area();


    let mut positions_op = path_renderers.iter_mut().map(|p| p.next()).collect::<Vec<_>>();
    let mut positions = positions_op.iter().map(|&p| p.unwrap()).collect::<Vec<_>>();

    while positions_op.iter().any(|p| p != &None) {
        vis_root.fill(&WHITE)?;
        let panels = vis_root.split_evenly((n_rows, n_cols));

        positions = positions.iter().zip(positions_op).map(|(&old, new)| new.unwrap_or(old) ).collect();

        // iteration over each path
        for (panel, pos) in panels.iter().zip(positions.iter()) {
            let mut graphic = ChartBuilder::on(&panel)
                .margin(1)
                .build_cartesian_2d(-1.1..1.1, -1.1..1.1 )?;

            graphic.draw_series( 
                embedding.iter()
                    .map(|(_, &v)| v)
                    .map(|p| Circle::new( 
                    (p[0], p[1]), 
                    2,
                    BLACK.filled()
                    )) 
            )?;

            // drawing edges
            for e in cplx.factor_pool.edge_iter()
                .map(|c| {
                    let r = &cplx.factor_pool.vertices_list();
                    let [x, y] = c.edge();
                    [embedding[&r[x]], embedding[&r[y]]]
                })
                .map(|[x, y]| [(x[0], x[1]), (y[0], y[1])] )
            {
                graphic.draw_series( LineSeries::new( 
                    e.into_iter(), 
                    BLACK.stroke_width(3)
                ))?;
            }

            // drawing robots
            graphic.draw_series(
                pos.iter().enumerate().map(|(i, p)| Circle::new( 
                    (p[0], p[1]), 
                    5,
                    VulcanoHSL::get_color(i as f64 / {N-1} as f64).filled()
                ))
            )?;
            graphic.draw_series(
                pos.iter().map(|p| Circle::new( 
                    (p[0], p[1]), 
                    5,
                    BLACK
                ))
            )?;
        }

        positions_op = path_renderers.iter_mut().map(|p| p.next()).collect::<Vec<_>>();

        vis_root.present()?;
    }

    for _ in 0..5 {
        vis_root.present()?;
    }

    Ok(())
}


pub fn draw_path_with_endpoints<const N: usize>(path: &cubical::EdgePath, start: &[usize], end: &[usize], embedding: &HashMap<SimpleVertex, ConstVector<f64, 2>>, graph: &SimpleGraph, graph_name: &str, bin_name: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut path_renderer = PathRenderer::<N>::new(
            path, 
            Some( start ),
            &embedding, 
            graph.vertices_list()
    );

    let into_coordinate = |cube_vertex: &[usize]| -> [ConstVector<f64, 2>; N] {
        let mut out = [ConstVector::<f64, 2>::zero(); N];
        for i in 0..N {
            out[i] = embedding[&graph.vertices_list()[cube_vertex[i]]]
        }
        out
    };
    let start = into_coordinate(start);
    let end = into_coordinate(end);

    // setting up the visualization tools
    let vis_root = BitMapBackend::gif(
        &format!("program_outputs/{bin_name}==={graph_name}==={N}_points.gif"), 
        (900, 300), 
        FLAME_DELAY
    )?.into_drawing_area();


    while let Some( position ) = path_renderer.next() {
        vis_root.fill(&WHITE)?;
        let panels = vis_root.split_evenly((1, 3));

        // drawing the three panels. 
        // the first one and the third one are the endpoints
        // the second one is the current position.
        for (panel, pos) in panels.iter().zip([start, position, end]) {
            let mut graphic = ChartBuilder::on(&panel)
                .margin(1)
                .build_cartesian_2d(-1.1..1.1, -1.1..1.1 )?;

            graphic.draw_series( 
                embedding.iter()
                    .map(|(_, &v)| v)
                    .map(|p| Circle::new( 
                    (p[0], p[1]), 
                    2,
                    BLACK.filled()
                    )) 
            )?;

            // drawing edges
            for e in graph.edge_iter()
                .map(|c| {
                    let r = &graph.vertices_list();
                    let [x, y] = c.edge();
                    [embedding[&r[x]], embedding[&r[y]]]
                })
                .map(|[x, y]| [(x[0], x[1]), (y[0], y[1])] )
            {
                graphic.draw_series( LineSeries::new( 
                    e.into_iter(), 
                    BLACK.stroke_width(3)
                ))?;
            }

            // drawing robots
            graphic.draw_series(
                pos.iter().enumerate().map(|(i, p)| Circle::new( 
                    (p[0], p[1]), 
                    5,
                    VulcanoHSL::get_color(i as f64 / {N-1} as f64).filled()
                ))
            )?;
            graphic.draw_series(
                pos.iter().map(|p| Circle::new( 
                    (p[0], p[1]), 
                    5,
                    BLACK
                ))
            )?;
        }

        vis_root.present()?;
    }

    for _ in 0..5 {
        vis_root.present()?;
    }

    Ok(())
}

pub fn draw_edge_path<const N: usize>(path: &CubicPath<N>, bin_name: &str, graph_name: &str, embedding: &[ConstVector<f64, 2>], graph: &RawSimpleGraph) -> Result<(), Box<dyn std::error::Error>> {
    // setting up the visualization tools
    let vis_root = BitMapBackend::gif(
        &format!("program_outputs/{bin_name}==={graph_name}==={N}_points.gif"), 
        (900, 300), 
        FLAME_DELAY
    )?.into_drawing_area();

    vis_root.fill(&WHITE)?;

    let start = path.start;

    let end = path.end;

    let mut pos = start;
    
    let panels = vis_root.split_evenly((1, 3));

    for (panel, pos) in panels.iter().zip([start, pos, end]) {
        let mut graphic = ChartBuilder::on(&panel)
            .margin(1)
            .build_cartesian_2d(-1.1..1.1, -1.1..1.1 )?;

        graphic.draw_series( 
            embedding.iter()
                .map(|p| Circle::new( 
                (p[0], p[1]), 
                2,
                BLACK.filled()
                )) 
        )?;

        // drawing edges
        for e in graph.edge_iter()
            .map(|c| {
                let [i, j] = c.edge();
                let (x, y) = (embedding[i], embedding[j]);
                [(x[0], x[1]), (y[0], y[1])]
            })
        {
            graphic.draw_series( LineSeries::new( 
                e.into_iter(), 
                BLACK.stroke_width(3)
            ))?;
        }

        // drawing robots
        graphic.draw_series(
            pos.iter().enumerate().map(|(i, &p)| Circle::new( 
                (embedding[p][0], embedding[p][1]), 
                5,
                VulcanoHSL::get_color(i as f64 / {N-1} as f64).filled()
            ))
        )?;
        graphic.draw_series(
            pos.iter().map(|&p| Circle::new( 
                (embedding[p][0], embedding[p][1]), 
                5,
                BLACK
            ))
        )?;
    }

    vis_root.present()?;

    for motion in path.iter() {
        panels[1].fill(&WHITE)?;

        motion.act_on(&mut pos);

        let mut graphic = ChartBuilder::on(&panels[1])
        .margin(1)
        .build_cartesian_2d(-1.1..1.1, -1.1..1.1 )?;

        graphic.draw_series( 
            embedding.iter()
                .map(|p| Circle::new( 
                (p[0], p[1]), 
                2,
                BLACK.filled()
                )) 
        )?;

        // drawing edges
        for e in graph.edge_iter()
            .map(|c| {
                let [i, j] = c.edge();
                let (x, y) = (embedding[i], embedding[j]);
                [(x[0], x[1]), (y[0], y[1])]
            })
        {
            graphic.draw_series( LineSeries::new( 
                e.into_iter(), 
                BLACK.stroke_width(3)
            ))?;
        }

        // drawing robots
        graphic.draw_series(
            pos.iter().enumerate().map(|(i, &p)| Circle::new( 
                (embedding[p][0], embedding[p][1]), 
                5,
                VulcanoHSL::get_color(i as f64 / {N-1} as f64).filled()
            ))
        )?;
        graphic.draw_series(
            pos.iter().map(|&p| Circle::new( 
                (embedding[p][0], embedding[p][1]), 
                5,
                BLACK
            ))
        )?;

        vis_root.present()?
    }

    for _ in 0..5 {
        vis_root.present()?;
    }

    Ok(())
}