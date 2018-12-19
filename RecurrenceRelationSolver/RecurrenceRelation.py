#!/usr/bin/env python3
# coding=utf-8
import logging
import re
import sympy
from sympy.solvers.solveset import linsolve

class RecurrenceVerificationFailed(Exception):
    """
    RecurrenceVerificationFailed will be thrown when verification fails with builtin
    sympy method.
    """
    def __init__(self, reason):
        """
        create RecurrenceVerificationFailed object
        Args:
            reason (string): Where the verification failed
        """
        self.reason = reason

class RecurrenceSolveFailed(Exception):
    """
    RecurrenceSolveFailed will be thrown when recurrence relation couldn't be solved fails
    """
    def __init__(self, reason):
        """
        create RecurrenceSolveFailed object
        Args:
            reason (string): Why the solve couldn't be performed
        """
        self.reason = reason

class RecurrenceRelation(object):
    """
    RecurrenceRelation object that contains a recurrence relations
    and initial conitions. Allows for solving and verifications
    """

    def __init__(self, recurrence, initialConditions):
        """
        create RecurrenceRelation object

        Args:
            recurrence (string): The recurrence as a string
            initialConditions (dict of string: string): The initial conditions
        """

        # contains the context for running sympy functions
        # on the recurrence relation
        self._sympy_context = {
            "s": sympy.Function("s"),
            "n": sympy.var("n", integer = True)
        }

        # Translate input string expression to sympy expression
        self._recurrence = self._to_sympy(recurrence)
        self._initialConditions = { k: self._to_sympy(v) for (k,v) in initialConditions.items() }
        
        # Solved values will be stored here in a bottom up dynamic programming manner
        self._solvedValues = dict(self._initialConditions)

        # Contains the closed from as calculated by our own algorithm
        self._closedForm = None

        self._logger = logging.getLogger(__name__)

    def _to_sympy(self, expr):
        """
        sympy represents powers not with ^ but with **.
        Translate from normal representation to sympy representation here
        
        Args:
            expr (string): string of an expression in normal format

        Returns:
            sympy expression: The  string parsed into a sympy expression
        """
        raw = re.sub("\^", "**", expr)

        return sympy.sympify(raw, self._sympy_context).expand()

    def _from_sympy(self, expr):
        """
        sympy represents powers not with ^ but with **.
        Translate from sympy representation to string representation here
        
        Args:
            expr (sympy expression): string of an expression in sympy format
        
        Returns:
            string: The sympy expression stringified
        """
        raw = str(expr)

        return re.sub("\*\*", "^", raw)

    def _getGeneralSolution(self, characteristicEq):
        """
        get the general solution given the characteristic equation.

        Args:
            characteristicEq (sympy expression): The characteristic equation
        
        Returns:
            sympy expression: The general solution
        """

        # get roots of characteristic equations and remove
        # the complex roots
        realRoots = { s:m for (s,m) in sympy.roots(characteristicEq).items() if sympy.I not in s.atoms() }
        # the sum of the multiplicity must be the same as the degree
        # else we can't solve the equation 

        if sum(realRoots.values()) != self._degree:
            msg = "The characteristic equation \"%s\" has the following real roots: %s, and the multiplicities is not the same as the degree" % (str(characteristicEq), str(realRoots))
            raise RecurrenceSolveFailed(msg)

        ctx = {
            "n": sympy.var("n", integer = True)
        }

        # Generate general solution
        generalSolutionTerms = []
        for i, (s,m) in enumerate(realRoots.items()):
            terms = []
            for j in range(0, m):
                varname = "p%d%d" % (i,j)
                ctx[varname] = sympy.var(varname)
                terms.append("(%s * n**%d)" % (varname, j))

            generalSolutionTerms.append("(%s)*(%s)**n" % ("+".join(terms), str(s)))


        return sympy.sympify('+'.join(generalSolutionTerms), ctx), ctx

    def _calculateClosedFromGeneralSolution(self, generalSolution, ctx):
        """
        get the closed form equation for a general solution.

        Args:
            generalSolution (sympy expression): The general solution
            ctx ((dict of string: sympy symbol):): The context of the general solution
        
        Returns:
            sympy expression: The closed form solved
        """

        # Create system of equations using initial conditions
        equations = []
        for i,c in self._initialConditions.items():
            eq = generalSolution - sympy.sympify("(%s)" % str(c))
            equations.append(eq.subs(ctx["n"], i))
        
        # Solve the system of equation
        solveSymbols = [ e for n, e in ctx.items() if n != "n" ]
        solutions = linsolve(equations, solveSymbols)
        
        if len(solutions) == 0:
            raise RecurrenceSolveFailed("No solution to the system of equations to find the alfas could be found.")

        solution = list(solutions)[0]

        # fill in the solution of the system
        solved = generalSolution
        for symbol, sub in zip(solveSymbols, list(solution)):
            solved = solved.subs(symbol, sub)

        return solved

    def _solveNonHomogeneous(self, generalSolution, ctx):
        """
        get the closed form equation for a non-homogeneous recurrence relation
        given the general solution for the associated homogeneous recurrence

        Args:
            generalSolution (sympy expression): The general solution for the associated homogeneous recurrence
            ctx ((dict of string: sympy symbol):): The context of the general solution
        
        Returns:
            sympy expression: The closed form solved
        """

        solveableRecurrence = self._recurrence - sympy.sympify("s(n)", self._sympy_context)
        # Guess an equation

        particularSolution = sympy.solve(solveableRecurrence)
        solution = particularSolution + generalSolution

        return self._calculateClosedFromGeneralSolution(solution, ctx)

    def _analyseExpression(self):
        """
        Analyse recurrence relation to determine certain properties
        
        Returns:
            tuple(int, sympy expr, sympy expr, bool): Contains information from the expression
                int: What degree the equation has
                sympy expr: the homogeneous part of the equation
                sympy expr: the non-homogeneous part of the equation
                bool: whether the recurrence is linear or not
        """

        homogenous = 0
        nonHomogenous = 0
        degree = 0
        linear = True

        s = self._sympy_context["s"]
        n = self._sympy_context["n"]
        i = sympy.Wild("i")

        for arg in self._recurrence.args:
            if arg.has(s):
                homogenous += arg
                rec_term = str(arg.find(s(n-i)))
                try:
                    term_degree = int(re.search(r'\d+', rec_term).group())
                    if term_degree > degree:
                        degree = term_degree
                    if linear and arg.has(sympy.Pow):
                        linear = False
                        break
                    if linear and arg.func == sympy.Mul:
                        s_args = arg.args
                        if (s_args[0].has(s) or s_args[0].has(s)) and (s_args[1].has(s) or s_args[1].has(n)):
                            linear = False
                            break
                except AttributeError:
                    raise RecurrenceSolveFailed("Failed parsing homogeneous term")
            else:
                nonHomogenous += arg

        return degree, homogenous, nonHomogenous, linear

    def _getDirectKey(self, expr):
        s = self._sympy_context["s"]
        # this is still way to hacky, maybe we have to check per step if we can do it or not...
        if expr.func == s:
            return expr.args[0].args[0]
        elif expr.args[1].func == s:
            return expr.args[1].args[0].args[0]


    def _dictBuilder (self, expr, dicts):
        s = self._sympy_context["s"]
        # if the function has an add, it means we can loop deeper
        if expr.func == sympy.Add :
            for arg in expr.args:
                self._dictBuilder(arg,dicts)
        elif expr.func == sympy.Mul:
            # we always get the value from a multiply function
            dicts[self._getDirectKey(expr)] =  expr.args[0]
        elif expr.func == s:
            # if we found the function directly, it means there is no extra mul, so we have a constant of 1
            dicts[self._getDirectKey(expr)] =  1
 
        return

    def _getCharacteristicEquation(self, expr):
        """
        Get the characteristic function with a given degree
        """    
        r = sympy.Symbol('r')

        #build the dictionary where we link every constant to every a_n value
        dicts = {}
        self._dictBuilder(self._recurrence, dicts)

        #start building our new expression!
        newExpr = r**self._degree
        for i in range(1, self._degree + 1):
            newExpr = newExpr - dicts.get(-i, 0) * r**(self._degree - i)
    
        return newExpr

    def _solve(self):
        """
        Solve the recurrence relation into a closed form

        Returns:
            String: The solved recurrence relation in string format
        """
        self._degree, homogenous, nonHomogenous, linear = self._analyseExpression()

        if not linear:
            raise RecurrenceSolveFailed("The equation is not linear")

        characteristicEq = self._getCharacteristicEquation(homogenous)
        generalSolution, ctx = self._getGeneralSolution(characteristicEq)

        solved = sympy.sympify("0")
        if nonHomogenous == 0:
            solved = self._calculateClosedFromGeneralSolution(generalSolution, ctx)
        else:
            solved = self._solveNonHomogeneous(generalSolution, ctx)

        return solved

    def solve(self):
        """
        Get the recurrence relation into a closed form

        Returns:
            String: The solved recurrence relation in string format
        """
        if self._closedForm is None:
            self._closedForm = self._solve()

        return self._from_sympy(self._closedForm)

    def _calculateValueFromRecurrence(self, n):
        """
        Get the nth value from the recurrence relation without solving it

        Args:
            n (int): The nth value to calculate
        
        Returns:
            sympy expression: The result

        """

        # Check if allready solved
        if n in self._solvedValues:
            return self._solvedValues[n]

        # Start solving from the next value that is not allready solved
        startSolving = max(self._solvedValues) + 1

        # remove the s(n) from the recurrence since we are not solving equal to 0
        # but equal to s(n)
        base = self._recurrence.subs(sympy.sympify("s(n)", self._sympy_context), 0)

        for i in range(startSolving, n + 1):
            eq = base
            # replace all function calls to itself with calculated values
            for j in range(1, self._degree + 1):
                replaceFunction = sympy.sympify("s(n-%d)" % j, self._sympy_context)
                replaceWith = self._solvedValues[i - j]
                eq = eq.subs(replaceFunction, replaceWith)

            # replace n with the current iteration, simplify the result and store it
            self._solvedValues[i] = eq.subs(self._sympy_context["n"], i).simplify()

        return self._solvedValues[n]


    def verify_range(self, upto, max_deviation):
        """
        Verify the solved recurrence relation with builtin sympy rsolve
        under the range [a,b]. If the difference for each value in
        the given range is larger then deviation an exception will
        be thrown giving the failed iteration.

        Args:
            upto (int): Upto what index to verify the solved formula with the recurrence
            max_deviation (float): What deviation to allow between the two solution and still be considered correct
        """

        # First solve the recurrence relation into closed form
        self.solve()

        print(self._closedForm)

        for i in range(min(self._initialConditions), upto + 1):
            iterative_result = self._calculateValueFromRecurrence(i).evalf()
            solved_result = self._closedForm.subs(self._sympy_context["n"], i).evalf()
            if abs(iterative_result - solved_result) >= max_deviation:
                msg = "Verification of solved recurrence failed at n = %d, our method gave: %f while iteration gave: %f." % (i, solved_result, iterative_result) 
                raise RecurrenceVerificationFailed(msg)
