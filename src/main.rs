extern crate ndarray;
extern crate rand;
extern crate rand_distr;

mod poisson_disk_sampling;

fn main() {
    // poisson_disk_sampling::recursive_neighbors(&vec![3, 3, 3, 3], 0);

    let mut indices = vec![3, 3, 3];
    // for outer_dim in 0..indices.len() {
    //     for _ in 0..3 {
    //         for inner_dim in outer_dim + 1..indices.len() {
    //             for _ in 0..3 {
    //                 for inner_inner_dim in inner_dim + 1..indices.len() {
    //                     for _ in 0..3 {
    //                         println!("{:?}", indices);
    //                         indices[inner_inner_dim] += 1;
    //                     }
    //                     indices[inner_inner_dim] -= 3;
    //                 }
    //                 indices[inner_dim] += 1;
    //             }
    //             indices[inner_dim] -= 3;
    //         }
    //         indices[outer_dim] += 1;
    //     }
    // }

    let mut new_indices;
    let one_dimension_count: usize = 3;
    let number_of_neighbors = one_dimension_count.pow(indices.len() as u32);
    for i in 0..number_of_neighbors {
        let mut index = i;
        new_indices = indices.clone();
        for dim in (1..indices.len()).rev() {
            let multiple = one_dimension_count.pow(dim as u32);
            let multiple_count = index / multiple;
            let dimension = indices.len() - 1 - dim;
            new_indices[dimension] = indices[dimension] + multiple_count;
            index -= multiple_count * multiple;
        }
        new_indices[indices.len() - 1] = indices[indices.len() - 1] + index % 3;
        println!("{:?}", new_indices);
    }

    // let points = poisson_disk_sampling::generate(vec![10.0, 10.0], 2.5, 30);

    // for point in points {
    //     println!("{}", point[0]);
    // }
}
