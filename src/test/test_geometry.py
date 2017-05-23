import numpy as np
import geometry as geo


class TestVector3D:
    vect = geo.Vector([1, 0, 0])

    def test_ndim(self):
        assert np.allclose(self.vect.vect, np.array([[1], [0], [0]]))

    def test_vect(self):
        assert self.vect.ndim == 3

    def test_call(self):
        assert np.allclose(self.vect(), np.array([1, 0, 0]))

    def test_length(self):
        assert self.vect.length == 1.0

    def test_translate(self):
        save = np.copy(self.vect())
        self.vect.translate([1, 1, 1])
        assert np.allclose(save, self.vect())

    def test_scale(self):
        vect = geo.Vector([1, -1, 2])
        vect.scale([2, 2, -2])
        assert np.allclose(np.array([2, -2, -4]), vect())

    def test_rotate(self):
        self.vect.rotate('y', 90)
        assert np.allclose(np.array([0, 0, -1]), self.vect())

    def test_scale_normal(self):
        vect = geo.NormalVector([1, 0, 0])
        vect.scale([2, 2, 2])
        assert np.allclose(np.array([1, 0, 0]), vect())

    def test_translate_point(self):
        point = geo.Point([1, -1, 2])
        point.translate([1, 1, 1])
        assert np.allclose(np.array([2, 0, 3]), point())


class TestScatterShape3D:
    def test_args(self):
        shape = geo.ScatterShape3D([0, 0, 0], [1, 1, 1])
        pt_min, pt_max = shape.args
        assert np.allclose(pt_min, np.array([0, 0, 0]))
        assert np.allclose(pt_max, np.array([1, 1, 1]))


class TestBrick:
    """A brick should basically always say True"""
    shape = geo.Brick([0, 0, 0], [1, 1, 1])

    def test_is_inside_single(self):
        assert self.shape.is_inside(np.array([.5, .5, .5]))

    def test_is_inside_array(self):
        assert all(self.shape.is_inside(np.array([[.5, .5, .5],
                                                  [.1, .1, .1]])))


class TestParallelepiped:
    """
            /|
        z  / |
        ^ /  |
        |/  /
        |  /
       /X /
      / |/_______ > x
     |  /
     | /
     |/

        O

    'X' is [0, 0, .5], right in the middle. Projection in the (normalized)
    parallelogram coordinates is [.5, .5, .5].
    'O' is [0, 0, -2], projected at [.5, .5, -2], hence outside the shape
    """
    shape = geo.Parallelepiped([-1, -1, -1],
                               [2, 0, 2],
                               [0, 2, 0],
                               [0, 0, 1])

    def test_mat(self):
        assert np.allclose(self.shape.mat, np.array([[2, 0, 0],
                                                     [0, 2, 0],
                                                     [2, 0, 1]]))

    def test_args(self):
        assert all(np.allclose(arg, targ)
                   for arg, targ in zip(self.shape.args,
                                        [[-1, -1, -1], [2, 0, 2],
                                         [0, 2, 0], [0, 0, 1]]))

    def test_bounding(self):
        assert np.allclose(self.shape.vects.pt_min(), np.array([-1, -1, -1]))
        assert np.allclose(self.shape.vects.pt_max(), np.array([1, 1, 2]))

    def test_project(self):
        assert np.allclose(
            self.shape.project([[0, 0, .5], [0, 0, -2]]),
            np.array([[.5, .5, .5], [.5, .5, -2]]).T)

    def test_is_inside(self):
        assert self.shape.is_inside([0, 0, .5])

    def test_is_not_inside(self):
        assert not self.shape.is_inside([0, 0, -2])
