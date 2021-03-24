/*
 * Copyright (C) 2016 (Salomon Sickert)
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

import acceptance.AcceptanceOmega;
import acceptance.AcceptanceType;
import automata.DA;
import automata.DASimplifyAcceptance;
import automata.HOAF2DA;
import com.google.common.collect.ImmutableList;
import de.tum.in.naturals.set.NatBitSet;
import jhoafparser.consumer.HOAIntermediateStoreAndManipulate;
import jhoafparser.transformations.ToStateAcceptance;
import owl.automaton.AutomatonUtil;
import owl.automaton.MutableAutomaton;
import owl.automaton.acceptance.GeneralizedRabinAcceptance;
import owl.automaton.acceptance.RabinAcceptance;
import owl.automaton.minimizations.GeneralizedRabinMinimizations;
import owl.automaton.minimizations.GenericMinimizations;
import owl.automaton.minimizations.MinimizationUtil;
import owl.automaton.transformations.RabinDegeneralization;
import owl.factories.Factories;
import owl.ltl.Formula;
import owl.ltl.LabelledFormula;
import owl.ltl.ROperator;
import owl.ltl.WOperator;
import owl.ltl.rewriter.RewriterFactory;
import owl.ltl.visitors.UnabbreviateVisitor;
import owl.run.env.DefaultEnvironment;
import owl.run.env.Environment;
import owl.translations.rabinizer.ImmutableRabinizerConfiguration;
import owl.translations.rabinizer.RabinizerBuilder;
import owl.translations.rabinizer.RabinizerState;
import parser.ast.Expression;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismSettings;

import java.util.Arrays;
import java.util.BitSet;
import java.util.EnumSet;

import static owl.automaton.output.HoaPrintable.Option.SIMPLE_TRANSITION_LABELS;

public class LTL2Rabinizer extends PrismComponent
{
	public LTL2Rabinizer(PrismComponent parent)
	{
		super(parent);
	}

	public DA<NatBitSet, ? extends AcceptanceOmega> rabinize(Expression ltl, AcceptanceType... allowedAcceptance) throws PrismException
	{
		getLog().println("\nBuilding Rabinizer automaton (for " + ltl + ")...");

		LabelledFormula ltlFormula = ExpressionToFormula.convert(ltl);
		LabelledFormula rewritten = RewriterFactory.apply(RewriterFactory.RewriterEnum.MODAL_ITERATIVE, ltlFormula);
		UnabbreviateVisitor visitor = new UnabbreviateVisitor(ROperator.class, WOperator.class);
		Formula formula = rewritten.accept(visitor);

		ImmutableRabinizerConfiguration configuration = ImmutableRabinizerConfiguration.builder().
				supportBasedRelevantFormulaAnalysis(true)
				.eager(true)
				.computeAcceptance(true).build();

		Environment environment = new DefaultEnvironment(false, false);
		Factories factories = environment.factorySupplier().getFactories(rewritten);

		MutableAutomaton<RabinizerState, GeneralizedRabinAcceptance> rabinizerAutomaton =
				RabinizerBuilder.rabinize(formula, factories, configuration, environment);
		MinimizationUtil.applyMinimization(rabinizerAutomaton, ImmutableList.of(GeneralizedRabinMinimizations::minimizeOverlap,
				GeneralizedRabinMinimizations::minimizeMergePairs,
				GenericMinimizations::removeTransientAcceptance,
				GeneralizedRabinMinimizations::minimizeComplementaryInf,
				GeneralizedRabinMinimizations::minimizeGloballyIrrelevant,
				GeneralizedRabinMinimizations::minimizeEdgeImplications,
				GeneralizedRabinMinimizations::minimizeSccIrrelevant,
				GeneralizedRabinMinimizations::minimizeTrivial,
				GeneralizedRabinMinimizations::minimizePairImplications,
				GeneralizedRabinMinimizations::minimizeMergePairs,
				GeneralizedRabinMinimizations::minimizeComplementaryInf,
				GeneralizedRabinMinimizations::minimizePairImplications,
				GeneralizedRabinMinimizations::minimizeEdgeImplications,
				GeneralizedRabinMinimizations::minimizeSccIrrelevant,
				GeneralizedRabinMinimizations::minimizeGloballyIrrelevant));

		boolean degeneralize = !Arrays.asList(allowedAcceptance).contains(AcceptanceType.GENERALIZED_RABIN);

		@SuppressWarnings("rawtypes")
		MutableAutomaton hoaResult;
		if (degeneralize) {
			getLog().println("Degeneralizing rabinizer automaton");
			MutableAutomaton<?, RabinAcceptance> degeneralized = RabinDegeneralization.degeneralize(rabinizerAutomaton);
			hoaResult = degeneralized;
			//noinspection unchecked
			AutomatonUtil.complete(hoaResult, Object::new, () -> {
				BitSet rejecting = new BitSet();
				degeneralized.getAcceptance().getPairs().forEach(pair -> {
					if (pair.hasFinite()) {
						rejecting.set(pair.getFiniteIndex());
					}
				});
				return rejecting;
			});
		} else {
			GeneralizedRabinAcceptance.normalize(rabinizerAutomaton);
			hoaResult = rabinizerAutomaton;
			//noinspection unchecked
			AutomatonUtil.complete(hoaResult, Object::new, () -> {
				BitSet rejecting = new BitSet();
				rabinizerAutomaton.getAcceptance().getPairs().forEach(pair -> rejecting.set(pair.getFiniteIndex()));
				return rejecting;
			});
		}

		HOAF2DA consumerDA = new HOAF2DA();
		hoaResult.toHoa(new HOAIntermediateStoreAndManipulate(consumerDA, new ToStateAcceptance()), EnumSet.of(SIMPLE_TRANSITION_LABELS));
		DA<NatBitSet, ? extends AcceptanceOmega> result = consumerDA.getDA();

		if (getSettings().getBoolean(PrismSettings.PRISM_NO_DA_SIMPLIFY)) {
			return result;
		}
		return DASimplifyAcceptance.simplifyAcceptance(this, result, allowedAcceptance);
	}
}
