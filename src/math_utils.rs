

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