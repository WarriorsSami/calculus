use calculus_lib::newton_secant;

fn f(x: f32) -> f32 {
    x * x.exp() - 1.0
}

fn main() {
    // newton method
    let mut x_prev = 1.0;
    let mut x_curr = 1.0;
    let eps = 1e-6;
    let mut i = 0;
    let it_max = 6;

    while f(x_curr).abs() > eps && i < it_max {
        let x_temp = x_curr;
        x_curr = newton_secant(x_prev, x_curr, f);
        x_prev = x_temp;
        i += 1;
    }

    match i <= it_max {
        true => println!("x = {:.6}", x_curr),
        false => println!("No solution found"),
    }
}
