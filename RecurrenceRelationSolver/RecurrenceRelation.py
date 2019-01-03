#!/usr/bin/env python3
# coding=utf-8
import logging
import re
import sympy
from sympy.solvers.solveset import linsolve

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

    def _to_sympy(self, expr):
        """
        sympy represents powers not with ^ but with **.
        Translate from normal representation to sympy representation here
        
        Args:
            expr (string): string of an expression in normal format

        Returns:
            sympy expression: The  string parsed into a sympy expression
        """
        raw = re.sub(r"\^", "**", expr)

        return sympy.sympify(raw, self._sympy_context).expand()

    def _from_sympy(self, expr):
        """
        sympy represents square roots as sqrt() but we require it to be ^(1/2)
        also powers are represented as ** by sympy but we need to have ^.
        This function takes a sympy expression and returns the appropiate string.
        
        Args:
            expr (sympy expression): string of an expression in sympy format
        
        Returns:
            string: The sympy expression stringified
        """

        expressionString = re.sub(r"\*\*", "^", str(expr))

        i = expressionString.find("sqrt")
        while i != -1:
            nestCount = 0
            endSqrtIndex = -1
            for j in range(i, len(expressionString)):
                if expressionString[j] == '(':
                    nestCount += 1
                elif expressionString[j] == ')':
                    if nestCount == 1:
                        endSqrtIndex = j + 1
                        break
                    nestCount -= 1

            sqrtExpr = expressionString[i+4:endSqrtIndex].strip()

            expressionString = "%s(%s^(1/2))%s" % (expressionString[0:i], sqrtExpr, expressionString[endSqrtIndex:])
            i = expressionString.find("sqrt")
        
        return expressionString

    def getRecurrence(self):
        """
        Get the recurrence

        Returns:
            string: The recurrence
        """
        return self._from_sympy(self._recurrence)

    def getLowerBoundDomain(self):
        """
        get the low bound where the recurrence is defined

        Returns:
            int: The start where the recurrence is defined
        """
        return min([int(k) for k,v in self._initialConditions.items()])

    def _set_free_variables_to_zero(self, solution):
        """
        Given a sympy solution with multiple solution fill out any free variables as 0.

        Args:
            solution (dict of sympy symbol: sympy expression): The solution for each symbol
                                                               where free variables have been
                                                               eliminated
        """

        # some symbols might not be in the solutions itself so add them here
        newSolution = { a:solution[a] if a in solution else a for sym, sol in solution.items() for a in sol.atoms(sympy.Symbol) }
        newSolution.update(solution)

        # Find all free variables and create a dict so we can easily replace them
        replaceDict = { sym:0 for sym, sol in newSolution.items() if sym == sol }

        return { symbol: expr.subs(replaceDict) for symbol, expr in newSolution.items() }


    def _getGeneralSolution(self, realRoots):
        """
        get the general solution given the roots of the characteristic equation.

        Args:
            realRoots (dict of sympy expr: int): The roots of the characteristic equation with multiplicities
        
        Returns:
            sympy expression: The general solution
        """

        ctx = {
            "n": self._sympy_context["n"]
        }

        # Generate general solution
        generalSolutionTerms = []
        for i, (s,m) in enumerate(realRoots.items()):
            terms = []
            for j in range(0, m):
                varname = "p_%d_%d" % (i,j)
                ctx[varname] = sympy.var(varname)
                terms.append("%s * n**%d" % (varname, j))

            generalSolutionTerms.append("(%s)*(%s)**n" % (" + ".join(terms), str(s)))


        return sympy.sympify(' + '.join(generalSolutionTerms), ctx), ctx

    def _calculateClosedFromGeneralSolution(self, generalSolution, ctx):
        """
        get the closed form equation for a general solution.

        Args:
            generalSolution (sympy expression): The general solution
            ctx (dict of string: sympy symbol): The context of the general solution
        
        Returns:
            sympy expression: The closed form solved
        """

        # Create system of equations using initial conditions
        equations = []
        for i,c in self._initialConditions.items():
            eq = generalSolution - sympy.sympify("(%s)" % str(c))
            equations.append(eq.subs(ctx["n"], i))

        logging.info("Solving the system of linear equations:")
        for e in equations:
            logging.info("\t%s" % str(e))
        
        # Solve the system of equation
        solve_symbols = [ e for n, e in ctx.items() if n != "n" ]
        solutions = linsolve(equations, solve_symbols)

        if len(solutions) == 0:
            raise RecurrenceSolveFailed("No solution to the system of equations to find the alfas could be found.")

        logging.info("Raw solutions with possibly free variables: %s" % str(solutions))

        # linsolve returns a set so we translate it to a dict for the _set_free_variables_to_zero function
        solution = { symbol:sol for symbol, sol in zip(solve_symbols, list(list(solutions)[0])) }

        logging.info("Dict solution with possibly free variables: %s" % str(solution))

        solution = self._set_free_variables_to_zero(solution)

        logging.info("Dict solution no free variables: %s" % str(solution))

        # fill in the solution of the system
        solved = generalSolution.subs(solution)

        return solved

    def _theorem6Classifier(self, expr, ctx, buckets):
        """
        Classify a single term in the non homogenous part of a recurrence relation.
        Classification means that that the term will be decomposed into the constant
        part, a part with a constant to the power of n and n to the power of a constant.
        

        Args:
            expr (sympy expression): The term to classify
            ctx (dict of string: sympy symbol): The context of the solution
            buckets (dict of sympy expression: dict of sympy expression: sympy symbol): Contains the classified terms
                                               in the form (the constant to the power n):
                                               (n to the power of a constant): (constant)
        """

        args = expr.args

        if expr.func != sympy.Mul:
            args = [expr]
        
        power = sympy.sympify("1", ctx)
        constant = sympy.sympify("1", ctx)
        poly = sympy.sympify("0", ctx)

        for a in args:
            if ctx["n"] not in a.atoms():
                constant = a
            elif a.func == sympy.Symbol:
                poly = sympy.sympify("1", ctx)
            elif a.func == sympy.Pow and a.args[0].func == sympy.Symbol:
                poly = a.args[1]
            elif a.func == sympy.Pow and a.args[1].func == sympy.Symbol:
                power = a.args[0]

        if power not in buckets:
            buckets[power] = {}

        buckets[power][poly] = constant

    def _theorem6SolutionBuilder(self, realroots, nonHomogenous, ctx):
        """
        Given the roots of the associated homogenous recurrence relation and the non homogenous part F(n)
        of the equation build a particular solution of the correct form.

        Args:
            realRoots (dict of sympy expr: int): The roots of the characteristic equation with multiplicities
            nonHomogenous (sympy expression): The part of the equation that makes the recurrence non homogenous
            ctx (dict of string: sympy symbol): The context of the solution to build
        
        Returns:
            sympy expression: The form of the particular solution
        """

        nonHomogenous = nonHomogenous.expand()

        buckets = {}
        if nonHomogenous.func != sympy.Add:
            self._theorem6Classifier(nonHomogenous, ctx, buckets)
        else:
            for arg in nonHomogenous.args:
                self._theorem6Classifier(arg, ctx, buckets)

        particularSolutionTerms = []
        for i, (power, polys) in enumerate(buckets.items()):
            multiplicity = realroots.get(power, 0)
            highestPoly = max(k for k, v in polys.items())
            terms = []
            for j in range(highestPoly, -1, -1):
                varname = "q_%d_%d" % (i,j)
                ctx[varname] = sympy.var(varname)
                terms.append("%s * n**%d" % (varname, j))

            particularSolutionTerms.append("n**%s*(%s)*(%s)**n" % (str(multiplicity), " + ".join(terms), str(power)))

        solutionOfCorrectForm = " + ".join(particularSolutionTerms)

        logging.info("Buckets from theorem6 classifier: %s" % str(buckets))
        logging.info("Solution must exist of form: %s" % solutionOfCorrectForm)

        return sympy.sympify(solutionOfCorrectForm, ctx)
    
    def _particularSolutionClassifier(self, expr, ctx, buckets):
        """
        Classify a single term of the recurrence relation where the particular solution has been substituted into 
        Classification means that that the term will be decomposed into buckets of x^n that contains a list of each
        term with the degree of that x^n.
        

        Args:
            expr (sympy expression): The term to classify
            ctx (dict of string: sympy symbol): The context of the solution
            buckets (dict of sympy expression: dict of sympy expression: sympy symbol): Contains the classified terms
        """
        args = expr.args

        if expr.func != sympy.Mul:
            args = [expr]
        
        matcherCtx = { "n": ctx["n"], "i": sympy.Wild("i") }
        matcher = sympy.sympify("n - i", matcherCtx)

        base = sympy.sympify("1", ctx)
        degree = sympy.sympify("0", ctx)
        term = expr

        for a in args:
            if a.func == sympy.Pow and len(a.args[0].atoms(sympy.Symbol)) == 0:
                base = a.args[0]
                m = a.args[1].match(matcher)
                degree = m[matcherCtx["i"]]
                term = term.subs(a, 1) 
                break

        if base not in buckets:
            buckets[base] = []

        buckets[base].append({
            "term": term,
            "degree": degree
        })

    def _solveParticularCoefficients(self, expr, ctx):
        """
        Solve for the coefficients of the particular solution.

        Args:
            expr (sympy expression): The recurrence with the particular solution substituted into it
            ctx (dict of string: sympy symbol): The context of the recurrence with the particular solution substituted into it

        Returns:
            (dict of sympy symbol: sympy expressiong): Dict containing values of the coefficients
        """
        expr = expr.expand()

        # built up dict containing
        buckets = {}
        if expr.func != sympy.Add:
            self._particularSolutionClassifier(expr, ctx, buckets)
        else:
            for arg in expr.args:
                self._particularSolutionClassifier(arg, ctx, buckets)

        logging.info("Buckets for solving coefficients of particular: ")
        for p,b in buckets.items():
            logging.info("\t%s: %s" % (str(p), str(b)))

        eqs = []
        for p,b in buckets.items():
            highestDegree = max([t["degree"] for t in b ])
            eq = 0
            for t in b:
                degreeDiff = highestDegree - t["degree"]
                eq = eq + p**degreeDiff * (t["term"])

            eqs.append(eq)

        logging.info("Each of the following equations has to be 0")
        for eq in eqs:
            logging.info("\t%s" % str(eq))

        # start solving that shit
        result = {}
        for eq in eqs:
            solveSymbols = [ s for s in eq.atoms(sympy.Symbol) if s != ctx["n"] ]
            res = sympy.solve_undetermined_coeffs(eq, solveSymbols, ctx["n"])
            if res is None:
                raise RecurrenceSolveFailed("Could not solve one of the sub equations of the filled in recurrence with the particular solution.")
                
            result.update(res)

        return result

    def _solveNonHomogeneous(self, realRoots, homogeneous, nonHomogenous, generalSolution, ctx):
        """
        get the closed form equation for a non-homogeneous recurrence relation
        given the general solution for the associated homogeneous recurrence

        Args:
            realRoots (dict of sympy expr: int): The roots of the characteristic equation with multiplicities
            homogeneous (sympy expression): The associated homogenous equation
            nonHomogenous (sympy expression): The part of the equation that makes the recurrence non homogenous
            generalSolution (sympy expression): The general solution for the associated homogeneous recurrence
            ctx (dict of string: sympy symbol): The context of the general solution
        
        Returns:
            sympy expression: The closed form solved
        """

        solveableRecurrence = self._recurrence - sympy.sympify("s(n)", self._sympy_context)

        logging.info("Brought all s(n) to one side for solving: %s" % str(solveableRecurrence))

        particularCtx = {
            "n": ctx["n"]
        }

        particularSolutionForm = self._theorem6SolutionBuilder(realRoots, nonHomogenous, particularCtx)

        logging.info("Particular solution is of form: %s" % str(particularSolutionForm))

        for i in range(0, self._degree + 1):
            particularReplaced = particularSolutionForm.subs(particularCtx["n"], particularCtx["n"] - i)
            replaceFunction = sympy.sympify("s(n-%d)" % i, self._sympy_context).simplify()
            solveableRecurrence = solveableRecurrence.subs(replaceFunction, particularReplaced)

        logging.info("Substituted the particular solution with the correct form into the recurrence: %s" % str(solveableRecurrence))

        # we replace every x^n with a positive constant so we can use some more sympy facilities.

        solution = self._solveParticularCoefficients(solveableRecurrence, particularCtx)

        logging.info("Solved for the variables: %s" % str(solution))

        # All solutions should not have any variables here any more if they have it means
        # the variable is free so we fill the free variables in with 0 here
        variables = self._set_free_variables_to_zero(solution)

        logging.info("Filled in any free variables as 0: %s" % str(variables))

        if solveableRecurrence.subs(variables).simplify() != 0:
            raise RecurrenceSolveFailed("Filling in solved coefficient for the recurrence where the particular solution has been fillled doesn't lead to a correct equation.")
 
        # Fill in the variables into the particular of the correct form
        particularSolution = particularSolutionForm.subs(variables)

        logging.info("Particular solution filled in with solved variables: %s" % str(particularSolution))

        generalSolution = (particularSolution + generalSolution).simplify()

        logging.info("Particular solution + general solution: %s" % str(generalSolution))

        solved = self._calculateClosedFromGeneralSolution(generalSolution, ctx)

        return solved

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
            if len(self._recurrence.args) == 1:
                arg = self._recurrence
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

        Returns:
            Sympy expression: The characteristic equation
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

        logging.info("Started solving recurrence relation: %s" % str(self._recurrence))

        self._degree, homogenous, nonHomogenous, linear = self._analyseExpression()

        msg = "homogenous" if nonHomogenous == 0 else "nonhomogenous"
        logging.info("Analyzation complete, It is a %s recurrence relation with degree %d" % (msg, self._degree))
        if nonHomogenous != 0:
            logging.info("The part that makes the recurrence nonhomogenous is: %s" % str(nonHomogenous))

        if not linear:
            raise RecurrenceSolveFailed("The equation is not linear")

        characteristicEq = self._getCharacteristicEquation(homogenous)
        logging.info("The characteristic equation is: %s" % str(characteristicEq))
       
        # get roots of characteristic equations and remove
        # the complex roots
        realRoots = { s:m for (s,m) in sympy.roots(characteristicEq).items() if sympy.I not in s.atoms() }
        logging.info("With roots: multiplicities: %s" % str(realRoots))

        # the sum of the multiplicity must be the same as the degree
        # else we can't solve the equation 
        if sum(realRoots.values()) != self._degree:
            msg = "The characteristic equation \"%s\" has the following real roots: %s, and the multiplicities is not the same as the degree" % (str(characteristicEq), str(realRoots))
            raise RecurrenceSolveFailed(msg)

        generalSolution, ctx = self._getGeneralSolution(realRoots)
        logging.info("The general solution has the form: %s" % str(generalSolution))

        solved = sympy.sympify("0")
        if nonHomogenous == 0:
            solved = self._calculateClosedFromGeneralSolution(generalSolution, ctx)
        else:
            solved = self._solveNonHomogeneous(realRoots, homogenous, nonHomogenous, generalSolution, ctx)
        
        logging.info("Solved raw: %s" % str(solved))
        solved = solved.simplify()
        logging.info("Solved simplified: %s" % str(solved))

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

    def calculateValueFromSolved(self, n):
        """
        Get the nth value from the solved recurrence relation

        Args:
            n (int): The nth value to calculate
        
        Returns:
            float: The result
        """
        self.solve()
        return self._closedForm.subs(self._sympy_context["n"], n).evalf(100)


    def calculateValueFromRecurrence(self, n):
        """
        Get the nth value from the recurrence relation without solving it

        Args:
            n (int): The nth value to calculate
        
        Returns:
            float: The result
        """

        # Check if allready solved
        if n in self._solvedValues:
            return self._solvedValues[n].evalf()

        # Start solving from the next value that is not allready solved
        startSolving = max(self._solvedValues) + 1

        for i in range(startSolving, n + 1):
            eq = self._recurrence
            # replace all function calls to itself with calculated values
            for j in range(1, self._degree + 1):
                replaceFunction = sympy.sympify("s(n-%d)" % j, self._sympy_context)
                replaceWith = self._solvedValues[i - j]
                eq = eq.subs(replaceFunction, replaceWith)

            # replace n with the current iteration, simplify the result and store it
            self._solvedValues[i] = eq.subs(self._sympy_context["n"], i).simplify()

        return self._solvedValues[n].evalf(100)

