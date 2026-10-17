"""Small invariant tests for leakage-sensitive code; no model performance tests."""
import unittest
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.base import clone
from workflow import LogCPM, FoldGeneFilter, make_pipeline, definitions, metrics, assert_group_split

class Invariants(unittest.TestCase):
    def test_cpm_uses_all_genes_and_is_sample_local(self):
        counts=np.array([[10,20,70],[2,8,10]],dtype=float)
        normalizer=LogCPM().fit(counts[:1])
        expected=np.log1p(counts/counts.sum(axis=1,keepdims=True)*1e6)
        np.testing.assert_allclose(normalizer.transform(counts),expected,rtol=1e-6)
        np.testing.assert_allclose(normalizer.transform(counts)[1],normalizer.transform(counts[1:])[0])
        with self.assertRaises(ValueError): normalizer.transform(np.zeros((1,3)))

    def test_validation_cannot_change_filter_or_scaler(self):
        train=np.array([[20,0,80],[50,0,50],[80,0,20],[30,0,70]],dtype=float)
        pipe=make_pipeline(definitions((0,1,2))['marker_logistic'][0])
        pipe.fit(train,np.array([0,0,1,1]))
        before=pipe.named_steps['filter'].indices_.copy()
        mean=pipe.named_steps['scale'].mean_.copy()
        self.assertNotIn(1,before)
        pipe.predict_proba(np.array([[1,1e8,1]],dtype=float))
        np.testing.assert_array_equal(before,pipe.named_steps['filter'].indices_)
        np.testing.assert_array_equal(mean,pipe.named_steps['scale'].mean_)
        clone(pipe)

    def test_grouping_with_repeated_patients(self):
        groups=np.repeat(np.arange(40),2)
        y=np.tile([0,0,1,1],20)
        for tr,va in StratifiedGroupKFold(5,shuffle=True,random_state=1).split(np.zeros((80,1)),y,groups):
            assert_group_split(tr,va,groups,y)

    def test_threshold_and_undefined_ppv(self):
        m=metrics(np.array([0,0,1,1]),np.array([0.1,0.2,0.2,0.3]))
        self.assertEqual(m['sensitivity'],0)
        self.assertEqual(m['specificity'],1)
        self.assertTrue(np.isnan(m['PPV']))
        self.assertEqual(metrics(np.array([0,1]),np.array([0.1,0.5]))['TP'],1)

if __name__=='__main__': unittest.main()
