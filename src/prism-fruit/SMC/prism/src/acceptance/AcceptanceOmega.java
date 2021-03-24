//==============================================================================
//
//	Copyright (c) 2014-
//	Authors:
//	* Joachim Klein <klein@tcs.inf.tu-dresden.de> (TU Dresden)
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

package acceptance;

import de.tum.in.naturals.set.NatBitSet;
import jdd.JDDVars;
import prism.PrismException;

import java.io.PrintStream;

/**
 * Generic interface for an omega-regular acceptance condition (NatBitSet-based).
 */
public interface AcceptanceOmega extends Cloneable
{
	/** Returns true if the bottom strongly connected component (BSSC)
	 *  given by bscc_states is accepting for this acceptance condition.
	 **/
	boolean isBSCCAccepting(NatBitSet bscc_states);

	/** Get the acceptance signature for state {@code stateIndex}.
	 **/
	String getSignatureForState(int stateIndex);

	/** Get the acceptance signature for state {@code stateIndex} (HOA format).
	 */
	String getSignatureForStateHOA(int stateIndex);

	/**
	 * Get a string describing the acceptance condition's size,
	 * i.e. "x Rabin pairs", etc.
	 */
	String getSizeStatistics();

	/** Returns the AcceptanceType of this acceptance condition */
	AcceptanceType getType();

	/** Returns the type of this acceptance condition as a String,
	 * i.e., "R" for Rabin.
	 * <br>
	 * Deprecated, use {@code getType().getNameAbbreviated()}
	 */
	@Deprecated String getTypeAbbreviated();

	/** Returns a full name for this acceptance condition
	 * <br>
	 * Deprecated, use {@code getType()} in String context or {@code getType().getName()}
	 */
	@Deprecated String getTypeName();

	/** Print the appropriate Acceptance (and potentially acc-name) header */
	void outputHOAHeader(PrintStream out);

	/** Make a copy of the acceptance condition. */
	AcceptanceOmega clone();

	/**
	 * Complement the acceptance condition if possible.
	 * @param numStates the number of states in the underlying model / automaton (needed for complementing NatBitSets)
	 * @param allowedAcceptance the allowed acceptance types that may be used for complementing
	 */
	AcceptanceOmega complement(int numStates, AcceptanceType... allowedAcceptance) throws PrismException;

	/** Abstract functor for use with the lift function. */
	abstract class LiftNatBitSet {
		public abstract NatBitSet lift(NatBitSet states);
	}

	/**
	 * Lift the state sets in the acceptance condition.
	 * For each state set {@code states} in the condition,
	 * {@code lifter.lift(states)} is called and the state set is
	 * replaced by the result.
	 **/
	void lift(LiftNatBitSet lifter);

	/**
	 * Convert this NatBitSet based acceptance condition to the corresponding BDD based acceptance condition.
	 * @param ddRowVars JDDVars of the row variables corresponding to the bits in the bitset
	 */
	AcceptanceOmegaDD toAcceptanceDD(JDDVars ddRowVars);

	/**
	 * Convert this acceptance condition to an AcceptanceGeneric condition.
	 */
	AcceptanceGeneric toAcceptanceGeneric();
}
