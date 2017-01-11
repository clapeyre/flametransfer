import pytest
from process_flames import *
 
# On passe de pytest.fixture() a pytest.yield_fixture()
@pytest.yield_fixture()
def simple_comme_bonjour():
    # tout ce qui est setup() va au dessus du yield. Ca peut etre vide.
    print('Avant !')
 
    # Ce qu'on yield sera le contenu du parametre. Ca peut etre None.
    yield ('pomme', 'banane')
 
    # Ce qu'il y a apres le yield est l'equivalent du tear down et peut etre
    # vide aussi
    print('Apres !')
 
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
