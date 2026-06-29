import numpy
import scipy

def line(
  x, 
  a
):
  
  return a * x

def gauss_left_tail(
  x, 
  A, 
  mu, 
  sigma, 
  k
):

  return numpy.where(
    x < mu,
    A * numpy.exp(k * (x - mu)),
    A * numpy.exp(-0.5 * ((x - mu) / sigma) ** 2),
  )

def crystal_ball_left(
  x, 
  A, 
  mu, 
  sigma, 
  alpha, 
  n
):
    """
      Standard Crystal Ball function (left-side power-law tail).
        x >= mu - alpha*sigma : Gaussian core
        x <  mu - alpha*sigma : power-law tail ~ (n/|alpha| - |alpha| - (x-mu)/sigma)^{-n}
      alpha > 0, n > 1
    """
    t = (x - mu) / sigma
    abs_a = numpy.abs(alpha)
    C = (n / abs_a) * (1.0 / (n - 1)) * numpy.exp(-0.5 * abs_a**2)
    D = numpy.sqrt(numpy.pi / 2) * (1 + scipy.special.erf(abs_a / numpy.sqrt(2)))
    tail_val = A * (n / abs_a)**n * numpy.exp(-0.5 * abs_a**2) * (n / abs_a - abs_a - t)**(-n)
    core_val = A * numpy.exp(-0.5 * t**2)

    return numpy.where(t < -abs_a, tail_val, core_val)

def crystal_ball_right(x, A, mu, sigma, alpha, n):
    """
    Crystal Ball with a right-side power-law tail.
    Reference: Laura++ manual, arXiv:1711.09854, eq. (108), sign-flipped.

    Gaussian core for t <= |alpha|, power-law tail for t > |alpha|,
    where t = (x - mu) / sigma.
    """
    t     = (x - mu) / sigma
    abs_a = numpy.abs(alpha)
    C     = (n / abs_a) ** n * numpy.exp(-0.5 * abs_a ** 2)
    B     = n / abs_a - abs_a
    return A * numpy.where(
        t <= abs_a,
        numpy.exp(-0.5 * t ** 2),      # Gaussian core
        C * (B + t) ** (-n),            # right power-law tail
    )


def crystal_ball_double(x, A, mu, sigma, alpha_l, n_l, alpha_r, n_r):
    """
    Double-sided Crystal Ball: power-law tails on both sides.
    Reference: arXiv:1711.09854 (Laura++), eq. (108) applied symmetrically.

    alpha_l, n_l : left  tail parameters
    alpha_r, n_r : right tail parameters
    """
    t      = (x - mu) / sigma
    abs_al = numpy.abs(alpha_l)
    abs_ar = numpy.abs(alpha_r)
    C_l    = (n_l / abs_al) ** n_l * numpy.exp(-0.5 * abs_al ** 2)
    B_l    = n_l / abs_al - abs_al
    C_r    = (n_r / abs_ar) ** n_r * numpy.exp(-0.5 * abs_ar ** 2)
    B_r    = n_r / abs_ar - abs_ar
    return A * numpy.where(
        t < -abs_al,
        C_l * (B_l - t) ** (-n_l),     # left  power-law tail
        numpy.where(
            t <= abs_ar,
            numpy.exp(-0.5 * t ** 2),  # Gaussian core
            C_r * (B_r + t) ** (-n_r), # right power-law tail
        ),
    )