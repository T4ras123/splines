use nannou::prelude::*;
use nannou::color::rgb_u32;

#[derive(Clone, Copy, Debug)]
struct Point {
    x: f32,
    y: f32,
}

impl Point {
    fn new(x: f32, y: f32) -> Self {
        Point { x, y }
    }
}

struct CubicSpline {

    points: Vec<Point>,
    y_values: Vec<f32>,
    lin_coeffs: Vec<f32>,
    quad_coeffs: Vec<f32>,
    cubic_coeffs: Vec<f32>,
}

impl CubicSpline {

    fn new(points: &[Point]) -> Self {

        if points.len() < 2 {
            panic!("Need at least 2 points to interpolate;")
        }
        let mut sorted_points = points.to_vec();
        sorted_points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap());

        let n = sorted_points.len();

        let y: Vec<f32> = sorted_points.iter().map(|p| p.y).collect();
        let x: Vec<f32> = sorted_points.iter().map(|p| p.x).collect();

        let y_values = y.clone();
        let mut linear_coefs = vec![0.0; n - 1];
        let mut cubic_coeffs = vec![0.0; n - 1];

        let mut h = vec![0.0; n - 1];
        for i in 0..n-1 {
            h[i] = x[i+1] - x[i];
            if h[i] == 0.0 {
                panic!("x values must be distinct");
            }
        }

        let mut alpha = vec![0.0; n-1];
        for i in 1..n-1 {
            alpha[i] = 3.0 * ((y_values[i + 1] - y_values[i]) / h[i] - (y_values[i] - y_values[i - 1]) / h[i - 1]);
        }


        let mut quadratic_coeffs = vec![0.0; n];
        let mut l = vec![0.0; n];
        let mut mu = vec![0.0; n];
        let mut z = vec![0.0; n]; 

        l[0] = 1.0;
        for i in 1..n - 1 {
            l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }
        
        l[n - 1] = 1.0;
        for j in (0..n - 1).rev() {
            quadratic_coeffs[j] = z[j] - mu[j] * quadratic_coeffs[j + 1];
            linear_coefs[j] = (y_values[j + 1] - y_values[j]) / h[j] - h[j] * (quadratic_coeffs[j + 1] + 2.0 * quadratic_coeffs[j]) / 3.0;
            cubic_coeffs[j] = (quadratic_coeffs[j + 1] - quadratic_coeffs[j]) / (3.0 * h[j]);
        }

        CubicSpline {
            points: sorted_points,
            y_values: y_values,
            lin_coeffs: linear_coefs,
            quad_coeffs: quadratic_coeffs,
            cubic_coeffs: cubic_coeffs,
        }
    }

    fn evaluate(&self, x: f32) -> f32 {
        let mut i = 0;
        while i < self.points.len() - 1 && x > self.points[i + 1].x {
            i += 1;
        }
        
        if i >= self.points.len() - 1 {
            return self.y_values[self.points.len() - 1];
        }
        if x < self.points[0].x {
            return self.y_values[0];
        }
        
        let dx = x - self.points[i].x;
        
        self.y_values[i] + dx * (self.lin_coeffs[i] + dx * (self.quad_coeffs[i] + dx * self.cubic_coeffs[i]))
    }
}


struct Model {

    control_points: Vec<Point>,
    spline: Option<CubicSpline>,
    dragging_point: Option<usize>,
    show_control_points: bool,
    resolution: usize,
}


fn model(app: &App) -> Model {
    app.new_window()
        .size(1600, 1200)
        .title("Cubic Spline Visualization")
        .view(view)
        .mouse_pressed(mouse_pressed)
        .mouse_released(mouse_released)
        .mouse_moved(mouse_moved)
        .key_pressed(key_pressed)
        .build()
        .unwrap();

    let control_points = vec![
        Point::new(-300.0, 0.0),
        Point::new(-150.0, 100.0),
        Point::new(0.0, -100.0),
        Point::new(150.0, 100.0),
        Point::new(300.0, 0.0),
    ];

    // Create the spline
    let spline = if control_points.len() >= 2 {
        Some(CubicSpline::new(&control_points))
    } else {
        None
    };

    Model {
        control_points,
        spline,
        dragging_point: None,
        show_control_points: true,
        resolution: 400, 
    }
}

fn update(_app: &App, model: &mut Model, _update: Update) {
    if model.control_points.len() >= 2 {
        model.spline = Some(CubicSpline::new(&model.control_points));
    } else {
        model.spline = None;
    }
}

fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    
    draw.background().color(rgb_u32(0x123456));

    if let Some(ref spline) = model.spline {
        let min_x = model.control_points.iter().map(|p| p.x).fold(f32::INFINITY, f32::min);
        let max_x = model.control_points.iter().map(|p| p.x).fold(f32::NEG_INFINITY, f32::max);
        
        let step = (max_x - min_x) / model.resolution as f32;
        let mut curve_points = Vec::with_capacity(model.resolution + 1);
        
        for i in 0..=model.resolution {
            let x = min_x + step * i as f32;
            let y = spline.evaluate(x);
            curve_points.push(pt2(x, y));
        }
        
        if curve_points.len() >= 2 {
            draw.polyline()
                .weight(3.0)
                .points(curve_points)
                .color(rgb_u32(0x00FFAA));
        }
    }

    if model.show_control_points {
        for (i, point) in model.control_points.iter().enumerate() {
            let is_selected = model.dragging_point == Some(i);
            let color = if is_selected { rgb_u32(0xFF3366) } else { rgb_u32(0xFFFFFF) };
            let size = if is_selected { 12.0 } else { 8.0 };
            
            draw.ellipse()
                .x_y(point.x, point.y)
                .radius(size)
                .color(color);
        }
    }

    let instructions = [
        "Click - Add Point",
        "Click+Drag - Move Point",
        "H - Toggle Control Points",
        "R - Reset Points",
        "C - Clear Points",
    ];

    for (i, text) in instructions.iter().enumerate() {
        draw.text(text)
            .x_y(-app.window_rect().w() / 2.0 + 100.0, app.window_rect().h() / 2.0 - 20.0 - i as f32 * 20.0)
            .color(WHITE)
            .font_size(14);
    }

    draw.to_frame(app, &frame).unwrap();
}

fn mouse_pressed(app: &App, model: &mut Model, button: MouseButton) {
    match button {
        MouseButton::Left => {
            let mouse_pos = app.mouse.position();
            let point = Point::new(mouse_pos.x, mouse_pos.y);
            
            let mut clicked_on_point = false;
            for (i, existing_point) in model.control_points.iter().enumerate() {
                let distance = ((existing_point.x - point.x).powi(2) + 
                                (existing_point.y - point.y).powi(2)).sqrt();
                if distance < 15.0 {
                    model.dragging_point = Some(i);
                    clicked_on_point = true;
                    break;
                }
            }
            
            if !clicked_on_point {
                model.control_points.push(point);
            }
        }
        _ => {}
    }
}

fn mouse_released(_app: &App, model: &mut Model, _button: MouseButton) {
    model.dragging_point = None;
}

fn mouse_moved(app: &App, model: &mut Model, pos: Vec2) {
    if let Some(idx) = model.dragging_point {
        model.control_points[idx].x = pos.x;
        model.control_points[idx].y = pos.y;
    }
}

fn key_pressed(app: &App, model: &mut Model, key: Key) {
    match key {
        Key::H => {
            model.show_control_points = !model.show_control_points;
        }
        Key::R => {
            model.control_points = vec![
                Point::new(-300.0, 0.0),
                Point::new(-150.0, 100.0),
                Point::new(0.0, -100.0),
                Point::new(150.0, 100.0),
                Point::new(300.0, 0.0),
            ];
            model.dragging_point = None;
        }
        Key::C => {
            model.control_points.clear();
            model.dragging_point = None;
        }
        Key::Escape => {
            app.quit();
        }
        _ => {}
    }
}

fn main() {
    nannou::app(model)
        .update(update)
        .run();
}
