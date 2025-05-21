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

#[derive(Clone, Copy, Debug, PartialEq)]
enum SplineType {
    Linear,
    Quadratic,
    Cubic,
}

struct Spline {
    points: Vec<Point>,
    spline_type: SplineType,
    a_coeffs: Vec<f32>,
    b_coeffs: Vec<f32>,
    c_coeffs: Vec<f32>,
    d_coeffs: Vec<f32>,
}

impl Spline {
    fn new(points: &[Point], spline_type: SplineType) -> Self {
        if points.len() < 2 {
            panic!("Need at least 2 points to interpolate;");
        }
        let mut sorted_points = points.to_vec();
        sorted_points.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap());

        let n = sorted_points.len();
        let x_coords: Vec<f32> = sorted_points.iter().map(|p| p.x).collect();
        let y_coords: Vec<f32> = sorted_points.iter().map(|p| p.y).collect();

        let mut h = vec![0.0; n - 1];
        for i in 0..n - 1 {
            h[i] = x_coords[i + 1] - x_coords[i];
            if h[i] == 0.0 {
                panic!("x values must be distinct for spline calculation.");
            }
        }

        let a_coeffs = y_coords.clone();
        let mut b_coeffs = vec![0.0; n - 1];
        let mut c_coeffs = vec![0.0; n - 1];
        let mut d_coeffs = vec![0.0; n - 1];

        match spline_type {
            SplineType::Linear => {
                for i in 0..n - 1 {
                    b_coeffs[i] = (a_coeffs[i + 1] - a_coeffs[i]) / h[i];
                }
            }
            SplineType::Quadratic => {
                if n >= 2 {
                    if h[0] != 0.0 {
                        b_coeffs[0] = (a_coeffs[1] - a_coeffs[0]) / h[0];
                        c_coeffs[0] = 0.0;
                    }

                    for i in 0..n - 2 {
                        let next_b;
                        if h[i] != 0.0 {
                            next_b = b_coeffs[i] + 2.0 * c_coeffs[i] * h[i];
                        } else {
                            next_b = 0.0;
                        }

                        if i + 1 < b_coeffs.len() {
                            b_coeffs[i + 1] = next_b;
                        }

                        if i + 1 < c_coeffs.len() && h[i + 1] != 0.0 {
                            c_coeffs[i + 1] = (a_coeffs[i + 2] - a_coeffs[i + 1] - b_coeffs[i + 1] * h[i + 1]) / (h[i + 1] * h[i + 1]);
                        }
                    }

                    if n == 2 && h[0] != 0.0 {
                        b_coeffs[0] = (a_coeffs[1] - a_coeffs[0]) / h[0];
                        c_coeffs[0] = 0.0;
                    }
                }
            }
            SplineType::Cubic => {
                let mut c_internal = vec![0.0; n];

                let mut alpha = vec![0.0; n - 1];
                for i in 1..n - 1 {
                    alpha[i] = 3.0 * ((a_coeffs[i + 1] - a_coeffs[i]) / h[i] - (a_coeffs[i] - a_coeffs[i - 1]) / h[i - 1]);
                }

                let mut l = vec![0.0; n];
                let mut mu = vec![0.0; n];
                let mut z = vec![0.0; n];

                l[0] = 1.0;

                for i in 1..n - 1 {
                    l[i] = 2.0 * (x_coords[i + 1] - x_coords[i - 1]) - h[i - 1] * mu[i - 1];
                    if l[i] == 0.0 {
                        panic!("Division by zero in cubic spline calculation (l[i])");
                    }
                    mu[i] = h[i] / l[i];
                    z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
                }

                l[n - 1] = 1.0;

                for j in (0..n - 1).rev() {
                    c_internal[j] = z[j] - mu[j] * c_internal[j + 1];

                    c_coeffs[j] = c_internal[j];
                    b_coeffs[j] = (a_coeffs[j + 1] - a_coeffs[j]) / h[j] - h[j] * (c_internal[j + 1] + 2.0 * c_internal[j]) / 3.0;
                    d_coeffs[j] = (c_internal[j + 1] - c_internal[j]) / (3.0 * h[j]);
                }
            }
        }

        Spline {
            points: sorted_points,
            spline_type,
            a_coeffs,
            b_coeffs,
            c_coeffs,
            d_coeffs,
        }
    }

    fn evaluate(&self, x: f32) -> f32 {
        if self.points.is_empty() {
            return 0.0;
        }
        if self.points.len() == 1 {
            return self.points[0].y;
        }

        let mut i = 0;
        while i < self.points.len() - 2 && x > self.points[i + 1].x {
            i += 1;
        }

        if x < self.points[0].x {
            if self.spline_type == SplineType::Linear && !self.b_coeffs.is_empty() {
                return self.a_coeffs[0] + self.b_coeffs[0] * (x - self.points[0].x);
            }
            return self.a_coeffs[0];
        }
        if x > self.points[self.points.len() - 1].x {
            let last_segment_idx = self.points.len().saturating_sub(2);
            if self.spline_type == SplineType::Linear && last_segment_idx < self.b_coeffs.len() {
                let dx_last = self.points[last_segment_idx + 1].x - self.points[last_segment_idx].x;
                if dx_last != 0.0 {
                    let slope = self.b_coeffs[last_segment_idx];
                    return self.a_coeffs[last_segment_idx + 1] + slope * (x - self.points[last_segment_idx + 1].x);
                }
            }
            return self.a_coeffs[self.points.len() - 1];
        }
        if x == self.points[self.points.len() - 1].x {
            return self.a_coeffs[self.points.len() - 1];
        }

        let dx = x - self.points[i].x;

        let mut val = self.a_coeffs[i];
        if i < self.b_coeffs.len() {
            val += self.b_coeffs[i] * dx;
        }

        if self.spline_type == SplineType::Quadratic || self.spline_type == SplineType::Cubic {
            if i < self.c_coeffs.len() {
                val += self.c_coeffs[i] * dx * dx;
            }
        }

        if self.spline_type == SplineType::Cubic {
            if i < self.d_coeffs.len() {
                val += self.d_coeffs[i] * dx * dx * dx;
            }
        }
        val
    }
}

struct Model {
    control_points: Vec<Point>,
    spline: Option<Spline>,
    dragging_point: Option<usize>,
    show_control_points: bool,
    resolution: usize,
    current_spline_type: SplineType,
}

fn model(app: &App) -> Model {
    app.new_window()
        .size(1600, 1200)
        .title("Spline Visualization")
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

    let current_spline_type = SplineType::Cubic;

    let spline = if control_points.len() >= 2 {
        Some(Spline::new(&control_points, current_spline_type))
    } else {
        None
    };

    Model {
        control_points,
        spline,
        dragging_point: None,
        show_control_points: true,
        resolution: 400,
        current_spline_type,
    }
}

fn update(_app: &App, model: &mut Model, _update: Update) {
    if model.control_points.len() >= 2 {
        model.spline = Some(Spline::new(&model.control_points, model.current_spline_type));
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

    let mut instructions = vec![
        "Click - Add Point",
        "Click+Drag - Move Point",
        "H - Toggle Control Points",
        "R - Reset Points",
        "C - Clear Points",
        "1 - Linear Spline",
        "2 - Quadratic Spline",
        "3 - Cubic Spline (Natural)",
    ];
    let current_spline_type_text = format!("Current Type: {:?}", model.current_spline_type);
    instructions.push(&current_spline_type_text);

    for (i, text) in instructions.iter().enumerate() {
        draw.text(text)
            .x_y(-app.window_rect().w() / 2.0 + 150.0, app.window_rect().h() / 2.0 - 20.0 - i as f32 * 20.0)
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
                let distance = ((existing_point.x - point.x).powi(2)
                    + (existing_point.y - point.y).powi(2))
                .sqrt();
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
        Key::Key1 => {
            model.current_spline_type = SplineType::Linear;
        }
        Key::Key2 => {
            model.current_spline_type = SplineType::Quadratic;
        }
        Key::Key3 => {
            model.current_spline_type = SplineType::Cubic;
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
