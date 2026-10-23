"""Independent rejection and arbitrary-center tests of coordinate certificates."""
import copy
import unittest

import sympy as sp

from verify_gram import verify_gram_solution


TETRA=[[0,1,2],[0,3,1],[0,2,3],[1,3,2]]
METRIC=[[1,sp.Rational(1,2),sp.Rational(1,2)],
        [sp.Rational(1,2),1,sp.Rational(1,2)],
        [sp.Rational(1,2),sp.Rational(1,2),1]]
POINTS=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]]


def tetra_payload(center):
    return {'metric':copy.deepcopy(METRIC),'coordinates':copy.deepcopy(POINTS),
            'center':list(center)}


class GramVerificationTests(unittest.TestCase):
    def assertCertified(self, report, location):
        self.assertTrue(report['verified'],report)
        self.assertTrue(report['certified'])
        self.assertEqual(report['status'],'CERTIFIED')
        self.assertEqual(report['center_location'],location)

    def assertRejected(self, report):
        self.assertFalse(report['verified'],report)
        self.assertEqual(report['status'],'UNVERIFIED')
        self.assertIn('reason',report)

    def test_interior_center_in_a_nonorthogonal_exact_basis(self):
        report=verify_gram_solution(TETRA,[[0,1,2,3]],[sp.Rational(3,8)],
                                   tetra_payload([sp.Rational(1,4)]*3))
        self.assertCertified(report,'interior')
        self.assertEqual(report['center_face_signs'],[1,1,1,1])

    def test_face_center_is_valid_for_a_sphere_cover(self):
        report=verify_gram_solution(TETRA,[[0,1,2],[3]],
                                   [sp.Rational(1,3),sp.Rational(2,3)],
                                   tetra_payload([sp.Rational(1,3),sp.Rational(1,3),0]))
        self.assertCertified(report,'boundary')
        self.assertEqual(report['center_face_signs'].count(0),1)

    def test_exterior_center_is_valid_for_a_sphere_cover(self):
        report=verify_gram_solution(TETRA,[[0],[1,2,3]],[6,3],
                                   tetra_payload([1,1,1]))
        self.assertCertified(report,'exterior')
        self.assertIn(-1,report['center_face_signs'])

    def test_center_classification_is_independent_of_face_orientation(self):
        report=verify_gram_solution([list(reversed(f)) for f in reversed(TETRA)],
                                   [[0],[1,2,3]],[6,3],tetra_payload([1,1,1]))
        self.assertCertified(report,'exterior')

    def test_claimed_group_radii_are_recomputed(self):
        for radii in ([6,4],[6,0],[-6,3]):
            self.assertRejected(verify_gram_solution(TETRA,[[0],[1,2,3]],radii,
                                                    tetra_payload([1,1,1])))

    def test_wrong_metric_or_duplicated_vertices_are_rejected(self):
        for key,value in [('metric',sp.diag(1,1,-1).tolist()),
                          ('coordinates',[[0,0,0],[0,0,0],[0,1,0],[0,0,1]])]:
            payload=tetra_payload([sp.Rational(1,4)]*3)
            payload[key]=value
            self.assertRejected(verify_gram_solution(TETRA,[[0,1,2,3]],
                                                    [sp.Rational(3,8)],payload))

    def test_floating_coefficients_and_malformed_payloads_are_rejected(self):
        payload=tetra_payload([sp.Rational(1,4)]*3)
        payload['metric'][0][0]=1.0
        for bad in (payload,{},None):
            self.assertRejected(verify_gram_solution(TETRA,[[0,1,2,3]],
                                                    [sp.Rational(3,8)],bad))

    def test_nonconvex_inward_cap_is_rejected_despite_exact_regular_faces(self):
        points=[[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                [0,0,1],[1,0,1],[1,1,1],[0,1,1],
                [sp.Rational(1,2),sp.Rational(1,2),1-sp.sqrt(2)/2]]
        faces=[[0,1,2,3],[0,4,5,1],[1,5,6,2],[2,6,7,3],[3,7,4,0],
               [4,7,8],[7,6,8],[6,5,8],[5,4,8]]
        center=[sp.Rational(1,2)]*3
        cap_radius=(sp.Rational(1,2)-sp.sqrt(2)/2)**2
        report=verify_gram_solution(faces,[list(range(8)),[8]],
                                   [sp.Rational(3,4),cap_radius],
                                   {'metric':sp.eye(3).tolist(),'coordinates':points,'center':center})
        self.assertRejected(report)
        self.assertIn('does not support',report['reason'])

    def test_primitive_embedding_is_checked(self):
        alpha=sp.Symbol('alpha')
        payload=tetra_payload([1,1,1])
        payload['metric'][0][0]=alpha**2/2
        report=verify_gram_solution(TETRA,[[0],[1,2,3]],[6,3],payload,
                                   primitive_symbol=alpha,primitive_polynomial=alpha**2-2,
                                   isolating_interval=(1,2))
        self.assertCertified(report,'exterior')
        self.assertRejected(verify_gram_solution(TETRA,[[0],[1,2,3]],[6,3],payload,
                                                 primitive_symbol=alpha,
                                                 primitive_polynomial=alpha**2-2,
                                                 isolating_interval=(2,3)))


if __name__=='__main__':unittest.main()
