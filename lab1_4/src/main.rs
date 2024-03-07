fn f(x: f32) -> f32 {
    (-x).exp()
}

fn main() {
    // newton method
    let mut x = 1.0;
    let eps = 1e-6;
    let mut i = 0;
    let it_max = 30;

    while (x - f(x)).abs() > eps && i < it_max {
        x = f(x);
        i += 1;
    }

    match i <= it_max {
        true => println!("x = {:.6}", x),
        false => println!("No solution found"),
    }
}
