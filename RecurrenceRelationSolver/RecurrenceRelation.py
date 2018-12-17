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

    def _getSolution(self, characteristicEq):
        """
        get the closed form equation given the characteristic equation.
        Get the fucking shit

        Args:
            characteristicEq (sympy expression): The characteristic equation
        
        Returns:
            sympy expression: The closed form solved
        """
        self._degree = 2
        self._initialConditions = {
            0: sympy.sympify(1),
            1: sympy.sympify(10)
        }
        characteristicEq = sympy.sympify("r**2 - 10*r + 1", { "r": sympy.var("r") })

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

        rawGeneralSolution = '+'.join(generalSolutionTerms)

        # Create system of equation using initial conditions
        equations = []
        for i,c in self._initialConditions.items():
            raw = rawGeneralSolution + ("-(%s)" % str(c))
            eq = sympy.sympify(raw, ctx).subs(ctx["n"], i)
            equations.append(eq)
        
        # Solve the system of equation
        solveSymbols = [ e for n, e in ctx.items() if n != "n" ]
        solutions = linsolve(equations, solveSymbols)
        
        if len(solutions) == 0:
            raise RecurrenceSolveFailed("No solution to the system of equations to find the alfas could be found.")

        solution = list(solutions)[0]

        # fill in the solution of the system
        solved = sympy.sympify(rawGeneralSolution, ctx)
        for symbol, sub in zip(solveSymbols, list(solution)):
            solved = solved.subs(symbol, sub)

        return solved

    def analyseExpression(self, expr):

        expandedTree = expr.expand()
        homogenous = 0
        nonHomogenous = 0
        degree = 0
        linear = True

        s = self._sympy_context["s"]
        n = self._sympy_context["n"]
        i = sympy.Wild("i")

        for arg in expandedTree.args:
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
                    print("Failed parsing homogeneous term")
            else:
                nonHomogenous += arg
        return degree, homogenous, nonHomogenous, linear


    def _solve(self):
        """
        Solve the recurrence relation into a closed form

        Returns:
            String: The solved recurrence relation in string format
        """
        self._degree, homogenous, nonHomogenous, linear = self.analyseExpression(self._recurrence)

        if not linear:
            msg = "The equation is not linear"
            raise RecurrenceSolveFailed(msg)

        return "Solved"

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

        # placeholder untill we actually get the degree but this should still work
        degree = len(self._initialConditions)

        # remove the s(n) from the recurrence since we are not solving equal to 0
        # but equal to s(n)
        base = self._recurrence.subs(sympy.sympify("s(n)", self._sympy_context), 0)

        for i in range(startSolving, n + 1):
            eq = base
            # replace all function calls to itself with calculated values
            for j in range(1, degree + 1):
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

        for i in range(min(self._initialConditions), upto + 1):
            iterative_result = self._calculateValueFromRecurrence(i).evalf()
            solved_result = iterative_result
            if abs(iterative_result - solved_result) >= max_deviation:
                msg = "Verification of solved recurrence failed at n = %d, our method gave: %f while iteration gave: %f." % (i, solved_result, iterative_result) 
                raise RecurrenceVerificationFailed(msg)
