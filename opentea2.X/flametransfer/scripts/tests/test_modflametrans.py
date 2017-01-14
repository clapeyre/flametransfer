from os.path import realpath, normpath
import pytest
from modflametrans import *
from XDRpy import init, finish, getDsout
 
# On passe de pytest.fixture() a pytest.yield_fixture()
@pytest.yield_fixture()
def dataset():
    # tout ce qui est setup() va au dessus du yield. Ca peut etre vide.
    init()
 
    # Ce qu'on yield sera le contenu du parametre. Ca peut etre None.
    yield getDsout()
 
    # Ce qu'il y a apres le yield est l'equivalent du tear down et peut etre
    # vide aussi
    finish()
 
def test_abs_path():
    from os import environ
    assert abs_path('..', 'TEST') == environ['C3SM_HOME']+'/library/flametransfer/scripts/TEST'

def test_proj_rel_path():
    assert proj_rel_path('GET', 'IN') == normpath(RUN_CURRENT + '/GET/IN')

def test_temp_rel_path():
    assert temp_rel_path('GET', 'IN') == normpath(TEMP + '/GET/IN')

def test_current_rel_path():
    assert current_rel_path('..', 'IN') == '../IN'

def test_current_common():
    assert current_common('file') == normpath('../'+COMMON+'/file')

#def test_get(simple_comme_bonjour):
#    element = get(simple_comme_bonjour, 0)
#    assert element == 'pomme'
# 
#def test_element_manquant(simple_comme_bonjour):
#    element = get(simple_comme_bonjour, 1000, 'Je laisse la main')
#    assert element == 'Je laisse la main'
# 
#def test_avec_echec(simple_comme_bonjour):
#    element = get(simple_comme_bonjour, 1000, 'Je laisse la main')
#    assert element == 'Je tres clair, Luc'
