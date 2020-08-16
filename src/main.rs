extern crate ndarray;
extern crate rand;
extern crate rand_distr;

mod poisson_disk_sampling;

fn main() {
    let points = poisson_disk_sampling::generate(&vec![10.0, 10.0], 2.5, 20);

    for point in points {
        println!("{:?}", point);
    }
}
