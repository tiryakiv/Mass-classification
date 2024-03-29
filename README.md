# Mass-classification

The codes for the paper entitled "Breast cancer mass classification using machine learning, binary-coded genetic algorithms, and an ensemble of deep transfer learning" are given in this repositorcy. The performance results were obtained by using mammograms from Breast Cancer Digital Repository (BCDR) [1-3].  01_DTL folder includes mass classification results on five-fold cross-validation data using ResNet50, NASNetLarge, and Xception deep transfer learning (DTL) methods. The performance results are for mammogram patch classification. 02_EDTL folder includes mass classification results on five-fold cross-validation data using ensembe of DTL (EDTL) methods. The performance results are for mammogram patch classification.

mass_class_edtl_valid_xval_perlesion_14feb23 includes mass lesion classification results using EDTL for lesions in the validation data.  

mass_class_edtl_valid_test_perlesion_15feb23 includes mass lesion classification results using EDTL for lesions in the test data. 

Please find the detailed information of our work from the following link. Also please cite our work if you find the content useful for your research:

https://academic.oup.com/comjnl/advance-article-abstract/doi/10.1093/comjnl/bxad046/7140287?utm_source=advanceaccess&utm_campaign=comjnl&utm_medium=email

Citation: Volkan Müjdat Tiryaki, Nedim Tutkun, Breast Cancer Mass Classification Using Machine Learning, Binary-Coded Genetic Algorithms and an Ensemble of Deep Transfer Learning, The Computer Journal, 2023;, bxad046, https://doi.org/10.1093/comjnl/bxad046


References:

1) 	(2012) , BCDR - Breast Cancer Digital Repository. [Online]. Available: http://bcdr.edu
2) 	Moura, D.C. et al. (2013) , Benchmarking Datasets for Breast Cancer Computer-Aided Diagnosis (CADx). , in Progress in Pattern Recognition, Image Analysis, Computer Vision, and Applications, pp. 326–333
3) 	Ramos-Pollán, R. et al. (2011) Discovering Mammography-based Machine Learning Classifiers for Breast Cancer Diagnosis. J. Med. Syst.,36, 2259–2269.
4)  L. Shen, L.R. Margoiles, J.H. Rothstein, E. Fluder, R. McBride, W. Sieh, Deep Learning to improve Breast cancer Detection on Screening Mammography, Sci. Rep. 9 (2019) 1–12. https://doi.org/https://doi.org/10.1038/s41598-019-48995-4.
5)  F. Chollet, Xception: Deep learning with depthwise separable convolutions, Proc. - 30th IEEE Conf. Comput. Vis. Pattern Recognition, CVPR 2017. 2017-Janua (2017) 1800–1807. https://doi.org/10.1109/CVPR.2017.195.
6)	K. He, X. Zhang, S. Ren, J. Sun, Deep residual learning for image recognition, in: Proc. IEEE Comput. Soc. Conf. Comput. Vis. Pattern Recognit., 2016: pp. 770–778. https://doi.org/10.1109/CVPR.2016.90.
7)	D.P. Kingma, J.L. Ba, Adam: A method for stochastic optimization, 3rd Int. Conf. Learn. Represent. ICLR 2015 - Conf. Track Proc. (2015) 1–15.
8)	F. Pedregosa, G. Varoquaux, A. Gramfort, V. Michel, B. Thirion, O. Grisel, M. Blondel, P. Prettenhofer, R. Weiss, V. Dubourg, J. Vanderplas, A. Passos, D. Cournapeau, M. Brucher, M. Perrot, E. Duchesnay, Scikit-learn: Machine Learning in Python, J. Mach. Learn. Res. 12 (2011) 2825–2830.
9)	B.W. Matthews, Comparison of the predicted and observed secondary structure of T4 phage lysozyme, Biochim. Biophys. Acta - Protein Struct. 405 (1975) 442–451. https://doi.org/10.1016/0005-2795(75)90109-9.
10) B. Zoph, V. Vasudevan, J. Shlens, Q. V. Le, Learning Transferable Architectures for Scalable Image Recognition, Proc. IEEE Comput. Soc. Conf. Comput. Vis. Pattern Recognit. (2018) 8697–8710. https://doi.org/10.1109/CVPR.2018.00907.
11) Jia Deng, Wei Dong, R. Socher, Li-Jia Li, Kai Li, Li Fei-Fei, ImageNet: A large-scale hierarchical image database, in: 2009 IEEE Conf. Comput. Vis. Pattern Recognit., IEEE, Miami, FL, USA, 2009: pp. 248–255. https://doi.org/10.1109/cvprw.2009.5206848.
