import re

def lex(monomial_1: tuple, monomial_2: tuple):
    for i in range(len(monomial_1)):
        if monomial_1[i] > monomial_2[i]:
            return True
        elif monomial_1[i] < monomial_2[i]:
            return False

    return True

def grlex(monomial_1: tuple, monomial_2: tuple):
    if sum(monomial_1) == sum(monomial_2):
        return lex(monomial_1, monomial_2)
    else:
        return sum(monomial_1) > sum(monomial_2)

def superscript(integer: str):
    substitution = integer.maketrans(''.join('0123456789'),
                                     ''.join('⁰¹²³⁴⁵⁶⁷⁸⁹'))
    
    return integer.translate(substitution)

def subscript(integer: str):
    substitution = integer.maketrans(''.join('0123456789'),
                                     ''.join('₀₁₂₃₄₅₆₇₈₉'))
    
    return integer.translate(substitution)

class Term:
    def __init__(self, coefficient: float = 0, variables: list[str] = [],
                 monomial: dict[str, int] = {}, order = lex, inplacer = None):
        for value in monomial.values():
            if value < 0:
                raise ValueError('Only positive exponents')

        if variables == []:
            variables = list(monomial.keys())

        if isinstance(inplacer, Term):
            coefficient = inplacer.coefficient
            variables = inplacer.variables
            monomial = inplacer.monomial
            order = inplacer.order

        self.coefficient = coefficient
        self.variables = variables
        self.monomial = monomial
        self.order = order
        self.clean()

    def __eq__(self, term):
        if not isinstance(term, Term):
            raise TypeError(str(type(term)) + 'is not comparable with Term')

        return self.monomial == term.monomial

    def __str__(self):
        if self.is_constant():
            return str(self['1'])

        string_representation = ''
        if self['1'] != 1:            
            string_representation = str(self['1'])
        if self['1'] == -1:
            string_representation = '-'

        for i in self.monomial:
            if self[i] == 1 :
                string_representation += i
            else:
                string_representation += i + superscript(str(self[i]))

        return string_representation

    def __len__(self):
        return len(self.monomial)
    
    def __getitem__(self, i: str):
        if i == '1':
            return self.coefficient

        return self.monomial[i]

    def __repr__(self):
        return self.monomial.__repr__()

    def __setitem__(self, i: str, newvalue: int):
        if i == '1':
            self.coefficient = newvalue

        self.monomial[i] = newvalue
        self.clean()

    def __gt__(self, term):
        if not isinstance(term, Term):
            raise TypeError(str(type(term)) + 'is not comparable with Term')

        if set(self.variables) != set(term.variables):
            raise ValueError('The terms must have the same variables')

        term_1 = [self[i] if i in self.monomial else 0 for i in self.variables]
        term_2 = [term[i] if i in term.monomial else 0 for i in self.variables]

        return self.order(term_1, term_2)

    def __neg__(self):
        return Term(coefficient = -self['1'], variables = self.variables,
                    monomial = self.monomial.copy(), order = self.order)

    def __mul__(self, term):
        if not isinstance(term, Term):
            raise TypeError(str(type(term)) + 'is not comparable with Term')

        new_monomial = self.monomial.copy()
        for i in term.monomial:
            if i in new_monomial:
                new_monomial[i] = self[i] + term[i]
            else:
                new_monomial[i] = term[i]

        return Term(coefficient = self['1'] * term['1'], order = self.order,
                    variables = self.variables, monomial = new_monomial)

    def __truediv__(self, term):
        if not isinstance(term, Term):
            raise TypeError(str(type(term)) + 'is not comparable with Term')

        new_monomial = self.monomial.copy()
        for i in term.monomial:
            if i in new_monomial:
                new_monomial[i] = self[i] - term[i]
            else:
                raise ValueError(str(self) + 'not divisible by' + str(term))

        return Term(coefficient = self['1'] / term['1'], order = self.order,
                    variables = self.variables, monomial = new_monomial)

    def is_constant(self):
        return len(self.monomial.keys()) == 0

    def from_string(string, order = lex):

        coefficient = re.findall('(?<=)[-\+]{0,1}[0-9]*', string)[0]

        expressions = re.findall('[a-zA-Z][_[0-9]*]{0,1}[\^[0-9]*]{0,1}', string)
        if coefficient == '' or coefficient == '+':
            coefficient = 1.0
        elif coefficient == '-':
            coefficient = -1.0
        else:
            coefficient = int(coefficient)

        monomial = {}
        for expression in expressions:
            exponent = re.findall('(?<=\^)[0-9]*', expression)
            variable = re.search('[a-zA-Z][_[0-9]*]{0,1}', expression)[0]
            if len(exponent) == 0:
                monomial[subscript(variable.replace('_', ''))] = 1
            else:
                monomial[subscript(variable.replace('_', ''))] = int(exponent[0])

        monomial = monomial
        variables = [re.findall('[a-zA-Z][_[0-9]*]{0,1}', expression)[0]
                              for expression in expressions]
        
        return Term(coefficient, variables, monomial, order)
            
        

    def divides(self, term):
        if not isinstance(term, Term):
            raise TypeError(str(type(term)) + 'is not comparable with Term')

        try:
            term / self
            return True
        except ValueError:
            return False

    def clean(self):
        if self['1'] == 0:
            self.monomial = {}

        monomial = self.monomial.copy()
        for i in self.monomial:
            if monomial[i] == 0:
                monomial.pop(i)

        self.monomial = monomial

class Polynomial:
    def __init__(self, variables: list[str], terms: list[Term] = [],
                 order = grlex, inplacer = None):
        if terms == []:
            terms = [Term(variables = variables, order = order)]

        if isinstance(inplacer, Term):
            terms = [inplacer]
            order = inplacer.order

        elif isinstance(inplacer, Polynomial):
            terms = inplacer.terms
            variables = inplacer.variables
            order = inplacer.order

        self.terms = terms
        self.variables = variables
        self.order = order
        self.clean()

    def clean(self):
        terms, new_terms = [], []
        [terms.append(term) for term in self.terms if term not in terms]
        for term in terms:
            sum_coefficient = sum([_.coefficient for _ in self.terms
                                   if _ == term])
            if sum_coefficient != 0:
                new_terms.append(Term(coefficient = sum_coefficient,
                                      monomial = term.monomial,
                                      variables = self.variables,
                                      order = self.order))

        self.terms = new_terms
        self.sort()

    @staticmethod
    def from_string(string, order = grlex):

        string = string.replace(" ", "")

        variables = list(set(re.findall('[a-zA-Z][_[0-9]*]{0,1}', string)))
        terms_str = re.split('([+-])', string)

        for i, term in enumerate(terms_str):

            if(term == '-' or term == '+'):
                terms_str[i] = terms_str[i] + terms_str[i + 1]
                del terms_str[ i + 1]

       # print(terms_str)

        terms = []
        for term in terms_str:
            terms.append(Term.from_string(term))
        variables = variables

        g = Polynomial(variables, terms, order)
        g.clean()

        return g

    def __eq__(self, poly):
        if isinstance(poly, Term):
           poly = Polynomial(variables = self.variables, inplacer = poly)

        if not isinstance(poly, Polynomial):
            raise TypeError(str(type(poly)) + 'is not comparable with Poynomial')

        return self.terms == poly.terms

    def __str__(self):
        if self.terms == []:
            return '0'

        string_representation = ''
        for term in self.terms:
            string_representation += str(term) + '+'

        return string_representation[0 : -1].replace('+-', '-')

    def __add__(self, poly):
        if isinstance(poly, Term):
            poly = Polynomial(variables = self.variables, inplacer = poly)

        if not isinstance(poly, Polynomial):
            raise TypeError(str(type(poly)) + 'is not comparable with Poynomial')

        h = Polynomial(terms = self.terms + poly.terms,
                       variables = self.variables, order = self.order)
        h.clean()
        return h

    def __neg__(self):
        return Polynomial(terms = [-term for term in self.terms],
                          variables = self.variables, order = self.order)

    def __sub__(self, poly):
        return self + (-poly)

    def __mul__(self, poly):
        if isinstance(poly, Term):
            poly = Polynomial(variables = self.variables, inplacer = poly)

        if not isinstance(poly, Polynomial):
            raise TypeError(str(type(poly)) + 'is not comparable with Poynomial')

        h = Polynomial(terms = [term1 * term2 for term1 in self.terms
                                for term2 in poly.terms],
                       variables = self.variables, order = self.order)
        h.clean()
        return h

    def __mod__(self, polys):
        _, r = GeneralizedDivAlg(self, polys)

        return r

    def sort(self):
        self.terms = [Term(coefficient = term.coefficient, order = self.order,
                           variables = self.variables, monomial = term.monomial)
                      for term in self.terms]
        self.terms = sorted(self.terms, reverse = True)

    def LT(self):
        self.sort()
        return Term(inplacer = self.terms[0])

    def LM(self):
        self.sort()
        term = Term(inplacer = self.terms[0])
        term.coefficient = 1
        return term


'''vars = ['x', 'y']

t1 = Term(2, (1, 3), vars)
t2 = Term(3, (0, 3), vars)
t3 = Term(1, (0, 1), vars)
t4 = Term(1, (1, 0), vars)

f = Polynomial([t1, t2, t3, t4], vars, order = grlex)
g = Polynomial([t1, t2, t3, t4], vars, order = grlex)

print(f)
print(g)
print(f + g)

t1 = Term(2, (1, 1), ['x', 'y'])
t2 = Term(3, (1, 0), ['x', 'y'])
print(t1 * t2)'''

def GeneralizedDivAlg(f: Polynomial, pols: list[Polynomial]):
    quotients = [Polynomial(variables = f.variables) for _ in range(len(pols))]
    r = Polynomial(variables = f.variables)
    p = Polynomial(variables = [], inplacer = f)

    while p != Polynomial(variables = f.variables):
        i = 0
        division_ocurred = False
        while i < len(pols) and division_ocurred == False:
            if pols[i].LT().divides(p.LT()):
                quotients[i] = quotients[i] + (p.LT() / pols[i].LT())
                p = p - (pols[i] * (p.LT() / pols[i].LT()))
                division_ocurred = True
            else:
                i = i + 1

        if division_ocurred == False:
            r = r + p.LT()
            p = p - p.LT()

    return quotients, r

def lcm(t1: Term, t2: Term):
    monomial = {}
    for variable in set(t1.variables + t2.variables):
        if variable in t1.monomial and variable in t2.monomial:
            monomial[variable] = max(t1[variable], t2[variable])
        elif variable in t1.monomial and variable not in t2.monomial:
            monomial[variable] = t1[variable]
        elif variable not in t1.monomial and variable in t2.monomial:
            monomial[variable] = t2[variable]

    return Term(coefficient = 1, monomial = monomial,
                variables = set(t1.variables + t2.variables), order = t1.order)

def S(f1: Polynomial, f2: Polynomial):
    gamma = lcm(f1.LM(), f2.LM())
    s1 = Polynomial(terms = [gamma / f1.LT()], variables = f1.variables,
                    order = f1.order)
    s2 = Polynomial(terms = [gamma / f2.LT()], variables = f2.variables,
                    order = f2.order)

    return s1 * f1 - s2 * f2

def GroebnerBasis(basis: list[Polynomial]):
    temp_basis = basis
    groebner_basis = basis
    for f_1 in temp_basis:
        for f_2 in temp_basis:
            residue = S(f_1, f_2) % temp_basis
            if residue != Polynomial(variables = f_1.variables):
                groebner_basis.append(residue)

    while temp_basis != groebner_basis:
        temp_basis = groebner_basis
        for f_1 in temp_basis:
            for f_2 in temp_basis:
                residue = S(f_1, f_2) % temp_basis
                if residue != Polynomial(variables = f_1.variables):
                    groebner_basis.append(residue)

    return groebner_basis

order = grlex

'''
f = Polynomial.from_string('x - y')
g = Polynomial.from_string('x + y - 1')

h = f * g

print(f * g)

q, r = GeneralizedDivAlg(h, [f])

print(q[0])'''

# BASES DE GROBNER

f1 = Polynomial.from_string("x + y", order = grlex)
print(f1 * f1)



