/*
 * Copyright (C) 2016  (Salomon Sickert)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package owl;

import it.unimi.dsi.fastutil.ints.IntAVLTreeSet;
import it.unimi.dsi.fastutil.ints.IntSortedSet;
import owl.ltl.BooleanConstant;
import owl.ltl.Conjunction;
import owl.ltl.Disjunction;
import owl.ltl.FOperator;
import owl.ltl.Formula;
import owl.ltl.GOperator;
import owl.ltl.LabelledFormula;
import owl.ltl.Literal;
import owl.ltl.ROperator;
import owl.ltl.UOperator;
import owl.ltl.WOperator;
import owl.ltl.XOperator;
import parser.ast.Expression;
import parser.ast.ExpressionBinaryOp;
import parser.ast.ExpressionLabel;
import parser.ast.ExpressionLiteral;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.type.TypeBool;
import prism.PrismLangException;

import java.util.ArrayList;
import java.util.List;

/*
 * Convert a property expression (an LTL formula) to classes understood by owl
 */
public final class ExpressionToFormula
{
	private ExpressionToFormula()
	{
	}

	public static LabelledFormula convert(Expression e)
	{
		if (e == null) {
			return null;
		}

		IntSortedSet variables = new IntAVLTreeSet();
		Formula formula = convert(e, variables);
		List<String> variableNames = new ArrayList<>();
		if (!variables.isEmpty()) {
			for (int i = 0; i <= variables.lastInt(); i++) {
				variableNames.add("L" + i);
			}
		}

		return LabelledFormula.create(formula, variableNames);
	}

	private static Formula convert(Expression e, IntSortedSet variables)
	{
		if (e == null) {
			return null;
		}

		if (e instanceof ExpressionTemporal) {
			return convertTemporal((ExpressionTemporal) e, variables);
		}
		if (e instanceof ExpressionBinaryOp) {
			return convertBinaryOp((ExpressionBinaryOp) e, variables);
		}
		if (e instanceof ExpressionUnaryOp) {
			return convertUnaryOp((ExpressionUnaryOp) e, variables);
		}
		if (e instanceof ExpressionLiteral) {
			return convertLiteral((ExpressionLiteral) e, variables);
		}
		if (e instanceof ExpressionLabel) {
			return convertLabel((ExpressionLabel) e, variables);
		}

		throw new IllegalArgumentException("Cannot convert expression: " + e.toString());
	}

	private static Formula convertTemporal(ExpressionTemporal e, IntSortedSet variables)
	{
		Formula ltl1 = convert(e.getOperand1(), variables);
		Formula ltl2 = convert(e.getOperand2(), variables);

		switch (e.getOperator()) {
		case ExpressionTemporal.P_X:
			return XOperator.create(ltl2);

		case ExpressionTemporal.P_U:
			return UOperator.create(ltl1, ltl2);

		case ExpressionTemporal.P_F:
			return FOperator.create(ltl2);

		case ExpressionTemporal.P_G:
			return GOperator.create(ltl2);

		case ExpressionTemporal.P_W:
			return WOperator.create(ltl1, ltl2);

		case ExpressionTemporal.P_R:
			return ROperator.create(ltl1, ltl2);

		default:
			throw new IllegalArgumentException("Cannot convert expression: " + e.toString());
		}
	}

	private static Formula convertBinaryOp(ExpressionBinaryOp e, IntSortedSet variables)
	{
		Formula ltl1 = convert(e.getOperand1(), variables);
		Formula ltl2 = convert(e.getOperand2(), variables);

		switch (e.getOperator()) {
		case ExpressionBinaryOp.IMPLIES:
			return Disjunction.create(ltl1.not(), ltl2);

		case ExpressionBinaryOp.IFF:
			return Conjunction.create(Disjunction.create(ltl1.not(), ltl2), Disjunction.create(ltl1, ltl2.not()));

		case ExpressionBinaryOp.OR:
			return Disjunction.create(ltl1, ltl2);

		case ExpressionBinaryOp.AND:
			return Conjunction.create(ltl1, ltl2);

		default:
			throw new IllegalArgumentException("Cannot convert expression: " + e.toString());
		}
	}

	private static Formula convertUnaryOp(ExpressionUnaryOp e, IntSortedSet variables)
	{
		Formula ltl1 = convert(e.getOperand(), variables);

		switch (e.getOperator()) {
		case ExpressionUnaryOp.NOT:
			return ltl1.not();

		case ExpressionUnaryOp.PARENTH:
			return ltl1;

		default:
			throw new IllegalArgumentException("Cannot convert expression: " + e.toString());
		}
	}

	private static Formula convertLiteral(ExpressionLiteral e, IntSortedSet variables)
	{
		if (!(e.getType() instanceof TypeBool)) {
			throw new IllegalArgumentException("Cannot convert expression: " + e.toString());
		}

		try {
			return BooleanConstant.get(e.evaluateBoolean());
		} catch (PrismLangException e1) {
			throw new IllegalArgumentException("Cannot convert expression: ", e1);
		}
	}

	private static Formula convertLabel(ExpressionLabel e, IntSortedSet variables)
	{
		int index = Integer.parseInt(e.getName().substring(1));
		variables.add(index);
		return new Literal(index, false);
	}
}
