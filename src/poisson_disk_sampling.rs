use ndarray::{Array, IxDyn};
use rand::prelude::*;
use rand_distr::{StandardNormal, Uniform};

macro_rules! clamp {
    ($value:expr, $min_value:expr, $max_value:expr) => {
        $value.max($min_value).min($max_value)
    };
}

type Point = Vec<f32>;

pub fn generate(extents: &Point, r: f32, k: i32) -> Vec<Point> {
    let dimensions = extents.len();
    let cell_size = r / (dimensions as f32).sqrt();
    let mut samples = Vec::<Point>::new();
    let mut active_list = Vec::<usize>::new();
    let background_grid_shape = point_add_scalar(&point_divide_scalar(&extents, cell_size), 1);
    let mut background_grid = Array::from_elem(background_grid_shape, -1 as isize);

    samples.push(random_point_in_square(&vec![0.0, 0.0], &extents));
    active_list.push(0);
    let indices = point_divide_scalar(&samples[0], cell_size);
    let ix_dyn = IxDyn(&indices);
    background_grid[ix_dyn] = 0;
    'outer: while !active_list.is_empty() {
        let active_point = &samples[active_list[active_list.len() - 1]];
        for _ in 0..k {
            let new_point = random_point_in_hollow_sphere(&active_point, r, r * 2.0);
            // TODO Pad the background grid instead?
            let clamped_point = clamp_point(&new_point, &vec![0.0, 0.0], &extents);
            if check_point_is_valid(&clamped_point, cell_size, &samples, &background_grid, r) {
                samples.push(clamped_point);
                active_list.push(samples.len() - 1);
                continue 'outer;
            }
        }
        // If no valid point could be generated
        active_list.pop();
    }

    samples
}

fn check_point_is_valid(
    point: &Vec<f32>,
    cell_size: f32,
    samples: &Vec<Vec<f32>>,
    background_grid: &Array<isize, IxDyn>,
    r: f32,
) -> bool {
    // These are the indices of the cell where the point would be located
    let mut indices = Vec::<usize>::with_capacity(point.len());
    for i in 0..point.len() {
        let index = (point[i] / cell_size).floor() as usize;
        indices.push(if index > 0 { index - 1 } else { index });
    }
    let indices = indices; // Remove mut

    let one_dimension_count: usize = 3;
    let number_of_neighbors = one_dimension_count.pow(indices.len() as u32);
    for i in 0..number_of_neighbors {
        let mut counter = i;
        let mut new_indices = indices.clone();
        for dim in (1..indices.len()).rev() {
            let multiple = one_dimension_count.pow(dim as u32);
            let multiple_count = counter / multiple;
            let dimension = indices.len() - 1 - dim;
            new_indices[dimension] = indices[dimension] + multiple_count;
            counter -= multiple_count * multiple;
        }
        new_indices[indices.len() - 1] = indices[indices.len() - 1] + counter % 3;

        let ix_dyn = IxDyn(&new_indices); // TODO add this to the ndarray doc
        let bg_index = background_grid[ix_dyn];
        if bg_index >= 0 {
            let background_point = &samples[bg_index as usize];
            if point_distance(point, background_point) < r / 2. {
                return false;
            }
        }
    }
    true
}

fn point_distance(p1: &Vec<f32>, p2: &Vec<f32>) -> f32 {
    assert!(p1.len() == p2.len());
    let dimensions = p1.len();
    let mut squared_sum = 0.0;
    for i in 0..dimensions {
        let difference = p1[i] - p2[i];
        squared_sum += difference * difference;
    }
    squared_sum.sqrt()
}

fn random_point_in_square(lower: &Vec<f32>, higher: &Vec<f32>) -> Vec<f32> {
    assert!(lower.len() == higher.len());
    let dimensions = lower.len();
    let mut point = Vec::<f32>::with_capacity(dimensions);
    for i in 0..dimensions {
        point.push(rand::random::<f32>() * (higher[i] - lower[i]) + lower[i]);
    }
    point
}

// https://en.wikipedia.org/wiki/N-sphere#Uniformly_at_random_within_the_n-ball
fn random_point_in_hollow_sphere(center: &Vec<f32>, r_min: f32, r_max: f32) -> Vec<f32> {
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
fn random_point_on_unit_n_sphere(dimensions: usize) -> Vec<f32> {
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

fn point_add_point(left: &Vec<f32>, right: &Vec<f32>) -> Vec<f32> {
    assert!(left.len() == right.len());
    let mut point = Vec::<f32>::with_capacity(left.len());
    for i in 0..left.len() {
        point.push(left[i] + right[i]);
    }
    point
}

fn point_add_scalar<Scalar>(left: &Vec<Scalar>, right: Scalar) -> Vec<Scalar>
where
    Scalar: std::ops::Add<Output = Scalar> + Copy,
{
    let mut point = Vec::<Scalar>::with_capacity(left.len());
    for i in 0..left.len() {
        point.push(left[i] + right);
    }
    point
}

fn point_mult_scalar(left: &Vec<f32>, right: f32) -> Vec<f32> {
    let mut point = Vec::<f32>::with_capacity(left.len());
    for i in 0..left.len() {
        point.push(left[i] * right);
    }
    point
}

fn point_divide_scalar(left: &Vec<f32>, right: f32) -> Vec<usize> {
    let mut point = Vec::<usize>::with_capacity(left.len());
    for i in 0..left.len() {
        point.push((left[i] / right).floor() as usize);
    }
    point
}
