use std::ops::{Add, Sub, Mul, Div};

#[derive(Copy, Clone)]
pub struct Conserved(pub f64, pub f64, pub f64);

#[derive(Copy, Clone)]
pub struct Primitive(pub f64, pub f64, pub f64);

#[derive(Copy, Clone)]
pub enum Direction { X, Y }




// ============================================================================
impl Direction
{
    fn dot(self, other: Direction) -> f64
    {
        match (self, other)
        {
            (Direction::X, Direction::X) => 1.0,
            (Direction::Y, Direction::Y) => 1.0,
            _ => 0.0,
        }
    }
}




// ============================================================================
impl Add<Primitive> for Primitive { type Output = Self; fn add(self, u: Primitive) -> Primitive { Primitive(self.0 + u.0, self.1 + u.1, self.2 + u.2) } }
impl Sub<Primitive> for Primitive { type Output = Self; fn sub(self, u: Primitive) -> Primitive { Primitive(self.0 - u.0, self.1 - u.1, self.2 - u.2) } }
impl Mul<f64> for Primitive { type Output = Primitive; fn mul(self, a: f64) -> Primitive { Primitive(self.0 * a, self.1 * a, self.2 * a) } }
impl Div<f64> for Primitive { type Output = Primitive; fn div(self, a: f64) -> Primitive { Primitive(self.0 / a, self.1 / a, self.2 / a) } }




// ============================================================================
impl Add<Conserved> for Conserved { type Output = Self; fn add(self, u: Conserved) -> Conserved { Conserved(self.0 + u.0, self.1 + u.1, self.2 + u.2) } }
impl Sub<Conserved> for Conserved { type Output = Self; fn sub(self, u: Conserved) -> Conserved { Conserved(self.0 - u.0, self.1 - u.1, self.2 - u.2) } }
impl Mul<f64> for Conserved { type Output = Conserved; fn mul(self, a: f64) -> Conserved { Conserved(self.0 * a, self.1 * a, self.2 * a) } }
impl Div<f64> for Conserved { type Output = Conserved; fn div(self, a: f64) -> Conserved { Conserved(self.0 / a, self.1 / a, self.2 / a) } }




// ============================================================================
impl Into<[f64; 3]> for Primitive { fn into(self) -> [f64; 3] { [self.0, self.1, self.2] } }
impl Into<[f64; 3]> for Conserved { fn into(self) -> [f64; 3] { [self.0, self.1, self.2] } }
impl From<[f64; 3]> for Primitive { fn from(a:  [f64; 3]) -> Primitive { Primitive(a[0], a[1], a[2]) } }
impl From<[f64; 3]> for Conserved { fn from(a:  [f64; 3]) -> Conserved { Conserved(a[0], a[1], a[2]) } }




// ============================================================================
impl Conserved
{
    pub fn density        (self)  -> f64 { self.0 }
    pub fn momentum_x     (self)  -> f64 { self.1 }
    pub fn momentum_y     (self)  -> f64 { self.2 }
    pub fn to_primitive   (self)  -> Primitive {
        Primitive(
            self.density(),
            self.momentum_x() / self.density(),
            self.momentum_y() / self.density())
    }
}




// ============================================================================
impl Primitive
{
    pub fn density   (self) -> f64 { self.0 }
    pub fn velocity_x(self) -> f64 { self.1 }
    pub fn velocity_y(self) -> f64 { self.2 }
    pub fn velocity  (self, direction: Direction) -> f64
    {
    	match direction {
    		Direction::X => self.velocity_x(),
    		Direction::Y => self.velocity_y(),
    	}
    }

    pub fn momentum_x(self) -> f64
    {
        self.density() * self.velocity_x()
    }

    pub fn momentum_y(self) -> f64
    {
        self.density() * self.velocity_y()
    }

    pub fn pressure(self, sound_speed_squared: f64) -> f64
    {
        self.density() * sound_speed_squared
    }

    pub fn to_conserved(self) -> Conserved
    {
        Conserved(
            self.density(),
            self.momentum_x(),
            self.momentum_y())
    }

    pub fn outer_wavespeeds(self, direction: Direction, sound_speed_squared: f64) -> (f64, f64)
    {
        let cs = sound_speed_squared.sqrt();
        let vn = self.velocity(direction);
        (vn - cs, vn + cs)
    }

    pub fn flux_vector(self, direction: Direction, sound_speed_squared: f64) -> Conserved
    {
        use Direction::*;
        let pg = self.pressure(sound_speed_squared);
        let vn = self.velocity(direction);
        let advective_term = self.to_conserved() * vn;
        let pressure_term = Conserved(0.0, pg * direction.dot(X), pg * direction.dot(Y));
        advective_term + pressure_term
    }
}




// ============================================================================
pub fn riemann_hlle(pl: Primitive, pr: Primitive, direction: Direction, sound_speed_squared: f64) -> Conserved {
    let ul = pl.to_conserved();
    let ur = pr.to_conserved();
    let fl = pl.flux_vector(direction, sound_speed_squared);
    let fr = pr.flux_vector(direction, sound_speed_squared);

    let (alm, alp) = pl.outer_wavespeeds(direction, sound_speed_squared);
    let (arm, arp) = pr.outer_wavespeeds(direction, sound_speed_squared);
    let ap = alp.max(arp).max(0.0);
    let am = alm.min(arm).min(0.0);

    (fl * ap - fr * am - (ul - ur) * ap * am) / (ap - am)
}
