use calculus_lib::Matrix;

fn main() {
    let mut A = Matrix::new(vec![
        vec![1.0, -1.0, 2.0],
        vec![2.0, 1.0, 1.0],
        vec![1.0, 2.0, -3.0],
    ]);

    for k in 1..A.n() {
        if A.at(k, k) != 0.0 {
            for i in k + 1..=A.n() {
                for j in k + 1..=A.m() {
                    A.set(i, j, A.at(i, j) - A.at(i, k) * A.at(k, j) / A.at(k, k));
                }
                A.set(i, k, 0.0);
            }
        } else {
            println!("Cannot triangulate matrix: zero on diagonal at ({}, {})", k, k);
            break;
        }
    }

    println!("{}", A);
}
