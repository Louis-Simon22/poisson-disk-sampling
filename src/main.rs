mod poisson_disk_sampling;

use poisson_disk_sampling::PoissonDiskSampling;

fn main() {
    let points = PoissonDiskSampling::generate_once(vec![10.0, 10.0], 2.5, 20, 2);

    for point in &points {
        println!("{:?}", point);
    }
}
