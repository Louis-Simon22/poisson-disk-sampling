extern crate rand;
extern crate rand_distr;
extern crate ndarray;

mod poisson_disk_sampling;

fn main() {
    let points = poisson_disk_sampling::generate(vec![10.0, 10.0], 2.5, 30);

    for point in points {
        println!("{}", point[0]);
    }
}
