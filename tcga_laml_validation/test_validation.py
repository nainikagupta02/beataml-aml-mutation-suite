"""Inference invariants only. No model training and no TCGA performance evaluation."""
import sys
sys.dont_write_bytecode=True
import unittest
from unittest.mock import patch
import joblib
import numpy as np
from validate import BEAT, OUT, predict_frozen, read_star, DATA, check_frozen_files

class FrozenInferenceTests(unittest.TestCase):
    def test_complete_gene_input_matches_original_pipeline(self):
        model=joblib.load(BEAT/'results/FLT3-ITD_elastic_net.joblib')
        width=model.named_steps['log_cpm'].n_features_in_
        X=np.ones((2,width),dtype=float)
        X[1,::3]=2
        p,z=predict_frozen(model,X,np.ones(width,dtype=bool))
        np.testing.assert_allclose(p,model.predict_proba(X)[:,1],atol=1e-12)

    def test_unavailable_genes_have_exact_zero_contribution_without_fit(self):
        model=joblib.load(BEAT/'results/NPM1_elastic_net.joblib')
        width=model.named_steps['log_cpm'].n_features_in_
        X=np.ones((2,width),dtype=float)
        selected=model.named_steps['filter'].indices_
        active=np.flatnonzero(model.named_steps['model'].coef_[0]!=0)
        missing_column=selected[active[0]]
        available=np.ones(width,dtype=bool); available[missing_column]=False
        X[:,missing_column]=0
        before=joblib.hash(model)
        with patch.object(type(model.named_steps['model']),'fit',side_effect=AssertionError('No fitting permitted')):
            p,z=predict_frozen(model,X,available)
        self.assertTrue((z[:,active[0]]==0).all())
        self.assertEqual(before,joblib.hash(model))
        expected=model.named_steps['model'].intercept_[0]+z@model.named_steps['model'].coef_[0]
        np.testing.assert_allclose(p,1/(1+np.exp(-expected)),rtol=1e-6)

    def test_star_column_is_integer_counts_and_ids_are_unique(self):
        d,values=read_star(DATA/'counts/0e7c1a03-0837-42b6-a60d-16be95603566.tsv')
        self.assertEqual(len(d),60660)
        self.assertTrue(d.stable_id.is_unique)
        self.assertTrue((values==np.floor(values)).all())
        self.assertIn('ENSG00000000003',set(d.stable_id))
        self.assertNotIn('N_unmapped',set(d.stable_id))

    def test_original_artifacts_are_unchanged(self):
        self.assertGreater(check_frozen_files(),0)

if __name__=='__main__': unittest.main()
