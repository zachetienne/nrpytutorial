from mpmath import mpf,mp
from UnitTesting.standard_constants import precision

# Dictionary of trusted values to be used throughout files.
# Standard precision and seed values are precision: 30, seed: 1234.
# Note that changing these may drastically change the calculated values.

mp.dps = precision
trusted_values_dict = dict()

# Generated on: 2019-06-22
# Reason for changing values: Brendan Drachler found a bug in IDValencia3velocityU[0].
trusted_values_dict['FishBoneMoncriefIDGlobals'] = {'hm1': mpf('-0.0429364995502448650557827344344'), 'rho_initial': mpf('2.26038704659529882516099442676'), 'IDalpha': mpf('0.646725235319608126854994121083'), 'IDgammaDD[0][0]': mpf('2.13479955812291560863615313951'), 'IDgammaDD[0][1]': mpf('-0.410524722407508344872095407149'), 'IDgammaDD[0][2]': mpf('0.801334226803907674342835380464'), 'IDgammaDD[1][0]': mpf('-0.410524722407508344872095407149'), 'IDgammaDD[1][1]': mpf('1.11858489740408392203774478605'), 'IDgammaDD[1][2]': mpf('-0.398998744578643500950609281897'), 'IDgammaDD[2][0]': mpf('0.801334226803907674342835380464'), 'IDgammaDD[2][1]': mpf('-0.398998744578643500950609281897'), 'IDgammaDD[2][2]': mpf('1.77285468301084672970131565121'), 'IDKDD[0][0]': mpf('-0.0168653943086464492703962463516'), 'IDKDD[0][1]': mpf('-0.0957514856135400691438745725662'), 'IDKDD[0][2]': mpf('-0.843969416553931128879391617796'), 'IDKDD[1][0]': mpf('-0.0957514856135400691438745725662'), 'IDKDD[1][1]': mpf('0.693548485453948604549775828158'), 'IDKDD[1][2]': mpf('-0.0868438208068101023930187731967'), 'IDKDD[2][0]': mpf('-0.843969416553931128879391617796'), 'IDKDD[2][1]': mpf('-0.0868438208068101023930187731967'), 'IDKDD[2][2]': mpf('-0.194285461315438191346692739218'), 'IDbetaU[0]': mpf('0.386431430889270419093131247622'), 'IDbetaU[1]': mpf('0.136618276576329210539871190106'), 'IDbetaU[2]': mpf('0.412837923505791089085073989069'), 'IDValencia3velocityU[0]': mpf('0.433292639899888785962776566577'), 'IDValencia3velocityU[1]': mpf('0.675772044179938876909747131769'), 'IDValencia3velocityU[2]': mpf('0.638351344526967177313628416414')}
