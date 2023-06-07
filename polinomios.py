from functools import cmp_to_key

class Term:

    def __init__(self, coefficient, monomyal_tuple, variables):

        self.coefficient = coefficient
        self.monomyal_tuple = monomyal_tuple
        self.variables = variables

    def __eq__(self, o):

        return self.monomyal_tuple == o.monomyal_tuple and self.coefficient == o.coefficient
    
    def __str__(self):
        
        if(self.is_constant()):
            return str(self.coefficient)
        
        string_representation = ""

        if(self.coefficient != 1):
            
            string_representation = str(self.coefficient)
            

        for i in range(len(self.monomyal_tuple)):

            if(self.monomyal_tuple[i] != 0):


                if(self.monomyal_tuple[i] == 1):
                    string_representation += self.variables[i]
                else:
                    string_representation += self.variables[i] + "^" + str(self.monomyal_tuple[i])

        return string_representation
    
    def is_constant(self):

        zero_tuple = tuple( [0 for _ in self.monomyal_tuple])

        return zero_tuple == self.monomyal_tuple

class Polynomial:

    def __init__(self, terms, variables, order):

        self.terms = terms
        self.variables = variables
        self.order = order

        self.clean()

    def clean(self): #HAY QUE PROBAR ESTE MÉTODO!!!

        monomyals = list(set([term.monomyal_tuple for term in self.terms]))

        new_terms = []

        for monomyal in monomyals:

            sum_coeff = 0

            for term in self.terms:

                if(term.monomyal_tuple == monomyal):

                    sum_coeff += term.coefficient

            if(sum_coeff != 0):
                new_terms.append(Term(sum_coeff, monomyal, self.variables))

        self.terms = new_terms

        self.sort()

    def __str__(self):

       # self.clean()

        string_representation = ""

        for term in self.terms:

            string_representation += str(term) + " + "

        return string_representation[0:-2]


    def __add__(self, o):

        h = Polynomial(self.terms + o.terms, self.variables, self.order)

        h.clean()

        return h
    
    @staticmethod
    def lex(monomyal_1, monomyal_2): #COMPARA monomyal_1 >= monomyal_2

        n = len(monomyal_1)

        for i in range(0, n):

            if(monomyal_1[i] > monomyal_2[i]):
                return True
            elif(monomyal_1[i] < monomyal_2[i]):
                return False

            
        return True
    
    @staticmethod
    def grlex(monomyal_1, monomyal_2):

        n = len(monomyal_1)

        sum_1 = 0
        sum_2 = 0

        for i in range(0, n):
            sum_1 += monomyal_1[i]
            sum_2 += monomyal_2[i]

        if(sum_1 == sum_2):
            return Polynomial.lex(monomyal_1, monomyal_2)
        else:
            return sum_1 > sum_2
        
    @staticmethod
    def compare(order):

        def inner(term_1, term_2):

            monomyal_1 = term_1.monomyal_tuple
            monomyal_2 = term_2.monomyal_tuple

            if monomyal_1 == monomyal_2:
                return 0
            elif order(monomyal_1, monomyal_2):
                return -1
            else:
                return 1
            
        return inner
    
    def sort(self):

        self.terms = sorted(self.terms, key=cmp_to_key(Polynomial.compare(self.order)))

vars = ["x", "y"]

t1 = Term(2, (1, 3), vars)
t2 = Term(3, (0, 3), vars)
t3 = Term(1, (0, 1), vars)
t4 = Term(1, (1, 0), vars)

f = Polynomial([t1, t2, t3, t4], vars, order = Polynomial.grlex)
g = Polynomial([t1, t2, t3, t4], vars, order = Polynomial.grlex)

print(f)
print(g)
print(f + g)


