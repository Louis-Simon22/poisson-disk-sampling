use ndarray::{Array, IxDyn};
use rand::prelude::*;
use rand_distr::{StandardNormal, Uniform};
use std::vec;

macro_rules! clamp {
    ($value:expr, $min_value:expr, $max_value:expr) => {
        $value.max($min_value).min($max_value)
    };
}

// TODO test this
// TODO extend the Point type with the add/divide etc methods as operations
// TODO Make this code generic using Index trait, implement one for vec

type Point = Vec<f32>;
pub struct PoissonDiskSampling {
    extents: Point,
    cell_size: f32,
    samples: Vec<Point>,
    background_grid: Array<Option<usize>, IxDyn>,
    active_list: Vec<usize>,
    minimum_distance: f32,
    maximum_iterations: u32,
    dimensions: u32,
}

impl PoissonDiskSampling {
    pub fn new(
        extents: Point,
        minimum_distance: f32,
        maximum_iterations: u32,
        dimensions: u32,
    ) -> PoissonDiskSampling {
        let cell_size = minimum_distance / (dimensions as f32).sqrt();
        let background_grid_shape =
            point_add_scalar(&point_integer_divide_scalar(&extents, cell_size), 1);
        let pds = PoissonDiskSampling {
            extents,
            cell_size,
            samples: Vec::<Point>::new(),
            background_grid: Array::from_elem(background_grid_shape, None),
            active_list: Vec::<usize>::new(),
            minimum_distance,
            maximum_iterations,
            dimensions,
        };
        pds
    }

    pub fn generate_once(
        extents: Point,
        minimum_distance: f32,
        maximum_iterations: u32,
        dimensions: u32,
    ) -> Vec<Point> {
        let mut pds =
            PoissonDiskSampling::new(extents, minimum_distance, maximum_iterations, dimensions);
        pds.generate();
        pds.samples
    }

    pub fn generate(&mut self) {
        self.samples
            .push(random_point_in_square(&vec![0.0, 0.0], &self.extents));
        self.active_list.push(0);
        let indices = point_integer_divide_scalar(&self.samples[0], self.cell_size);
        let ix_dyn = IxDyn(&indices);
        self.background_grid[ix_dyn] = Some(0);
        'outer: while !self.active_list.is_empty() {
            let active_point = &self.samples[self.active_list[self.active_list.len() - 1]];
            for _ in 0..self.maximum_iterations {
                let new_point = random_point_in_hollow_sphere(
                    &active_point,
                    self.minimum_distance,
                    self.minimum_distance * 2.0,
                );
                // TODO For performance maybe add padding to the background grid instead?
                let clamped_point = clamp_point(&new_point, &vec![0.0, 0.0], &self.extents);
                if self.check_point_is_valid(&clamped_point) {
                    self.samples.push(clamped_point);
                    self.active_list.push(self.samples.len() - 1);
                    continue 'outer;
                }
            }
            // If no valid point could be generated
            self.active_list.pop();
        }
    }

    fn check_point_is_valid(&self, point: &Point) -> bool {
        // These are the indices of the cell where the point would be located
        let mut cell_indices = Vec::<usize>::with_capacity(point.len());
        for i in 0..point.len() {
            let index = (point[i] / self.cell_size).floor() as usize;
            cell_indices.push(if index > 0 { index - 1 } else { index });
        }
        let cell_indices = cell_indices;

        // indices = [3, 5, 10]
        // dimension_basis = [9, 3, 1]
        // neighbor_count = pow(3, len(indices))
        // print(neighbor_count)
        // for i in range(neighbor_count):
        //     current_count = i
        //     new_indices = indices.copy()
        //     for dim,basis in enumerate(dimension_basis):
        //         dim_count = current_count // basis
        //         current_count -= dim_count * basis
        //         new_indices[dim] += dim_count
        //     print(new_indices)

        // TODO comment this code
        const ONE_DIMENSION_NEIGHBOR_CELLS_COUNT: usize = 3;
        let neighbor_cells_count = ONE_DIMENSION_NEIGHBOR_CELLS_COUNT.pow(self.dimensions);
        for i in 0..neighbor_cells_count {
            let mut counter = i;
            let mut new_indices = cell_indices.clone();
            for dim in (0..cell_indices.len()).rev() {
                let multiple = ONE_DIMENSION_NEIGHBOR_CELLS_COUNT.pow(dim as u32);
                let multiple_count = counter / multiple;
                let new_index = new_indices[dim] + multiple_count;
                new_indices[dim] = match new_index.checked_sub(1) {
                    Some(new_index) => new_index,
                    None => new_index
                };
                counter -= multiple_count * multiple;
            }
            // TODO new_indices is sometimes out of range (clamp it? Make the array bigger?)

            let ix_dyn = IxDyn(&new_indices); // TODO add this to the ndarray doc
            let bg_index_opt = self.background_grid[ix_dyn];
            match bg_index_opt {
                Some(bg_index) => {
                    let background_point = &self.samples[bg_index];
                    if point_distance(point, background_point) < self.minimum_distance / 2.0 {
                        return false;
                    }
                },
                None => {}
            }
        }
        true
    }
}

fn point_distance(p1: &Point, p2: &Point) -> f32 {
    assert!(p1.len() == p2.len());
    let dimensions = p1.len();
    let mut squared_sum = 0.0;
    for i in 0..dimensions {
        let difference = p1[i] - p2[i];
        squared_sum += difference * difference;
    }
    squared_sum.sqrt()
}

fn random_point_in_square(lower: &Point, higher: &Point) -> Point {
    assert!(lower.len() == higher.len());
    let dimensions = lower.len();
    let mut point = Vec::<f32>::with_capacity(dimensions);
    for i in 0..dimensions {
        point.push(rand::random::<f32>() * (higher[i] - lower[i]) + lower[i]);
    }
    point
}

// https://en.wikipedia.org/wiki/N-sphere#Uniformly_at_random_within_the_n-ball
fn random_point_in_hollow_sphere(center: &Point, r_min: f32, r_max: f32) -> Point {
    let dimensions = center.len();
    let point_on_unit_n_sphere = random_point_on_unit_n_sphere(dimensions);
    // TODO check that the radius of the unit sphere is correct
    //let unit_radius =
    let radius = rand::thread_rng().sample(Uniform::from(r_min..r_max));
    let point_in_n_sphere = point_mult_scalar(
        &point_on_unit_n_sphere,
        radius.powf(1.0 / dimensions as f32),
    );
    point_add_point(&point_in_n_sphere, center)
}

// https://en.wikipedia.org/wiki/N-sphere#Generating_random_points
fn random_point_on_unit_n_sphere(dimensions: usize) -> Point {
    let mut squared_sum: f32 = 0.0;
    let mut point = Vec::<f32>::with_capacity(dimensions);
    for _ in 0..dimensions {
        let normal_deviate = rand::thread_rng().sample(StandardNormal);
        point.push(normal_deviate);
        squared_sum += normal_deviate * normal_deviate;
    }
    let radius = squared_sum.sqrt();
    for i in 0..point.len() {
        point[i] /= radius;
    }
    point
}

fn clamp_point(point: &Point, lower: &Point, higher: &Point) -> Point {
    let mut clamped_point = Point::with_capacity(point.len());
    for i in 0..point.len() {
        clamped_point.push(clamp!(point[i], lower[i], higher[i]));
    }
    clamped_point
}

fn point_add_point(left: &Point, right: &Point) -> Point {
    assert!(left.len() == right.len());
    let mut point = Vec::<f32>::with_capacity(left.len());
    for i in 0..left.len() {
        point.push(left[i] + right[i]);
    }
    point
}

fn point_add_scalar<Scalar>(left: &Vec<Scalar>, right: Scalar) -> Vec<Scalar>
where
    Scalar: std::ops::Add<Output = Scalar> + Copy + std::fmt::Debug,
{
    let mut point = Vec::<Scalar>::with_capacity(left.len());
    for i in 0..left.len() {
        point.push(left[i] + right);
    }
    point
}

fn point_mult_scalar(left: &Point, right: f32) -> Point {
    let mut point = left.clone();
    for i in 0..left.len() {
        point[i] = left[i] * right;
    }
    point
}

fn point_integer_divide_scalar(left: &Point, right: f32) -> Vec<usize> {
    let mut point = Vec::<usize>::with_capacity(left.len());
    for i in 0..left.len() {
        point.push((left[i] / right) as usize);
    }
    point
}
