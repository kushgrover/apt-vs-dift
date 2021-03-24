package common;

import de.tum.in.naturals.Indices;
import explicit.MEC;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMaps;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntIterators;
import it.unimi.dsi.fastutil.ints.IntSet;
import prism.PrismUtils;

import java.util.Iterator;
import java.util.function.IntConsumer;
import java.util.function.IntUnaryOperator;

import static com.google.common.base.Preconditions.checkArgument;

public class FastUtils
{
	/**
	 * Feed all elements of {@code ints} to the {@code consumer} or every integer between {@code 0} inclusive and {@code length} exclusive, if {@code ints}
	 * is {@code null}.
	 */
	public static void forEach(IntIterable ints, int length, IntConsumer consumer)
	{
		if (ints == null) {
			for (int i = 0; i < length; i++) {
				consumer.accept(i);
			}
		} else {
			ints.forEach(consumer);
		}
	}

	public static boolean doublesAreClose(Int2DoubleMap d1, Int2DoubleMap d2, double precision, boolean absolute)
	{
		checkArgument(d1.size() == d2.size());

		Iterator<Int2DoubleMap.Entry> d1iterator = Int2DoubleMaps.fastIterator(d1);
		Iterator<Int2DoubleMap.Entry> d2Iterator = Int2DoubleMaps.fastIterator(d2);

		while (d1iterator.hasNext()) {
			assert d2Iterator.hasNext();
			Int2DoubleMap.Entry d1next = d1iterator.next();
			Int2DoubleMap.Entry d2next = d2Iterator.next();
			if (d1next.getIntKey() != d2next.getIntKey()) {
				return false;
			}
			if (!PrismUtils.doublesAreClose(d1next.getDoubleValue(), d2next.getDoubleValue(), precision, absolute)) {
				return false;
			}
		}
		return true;
	}

	public static IntUnaryOperator elementToIndexMap(IntSet set)
	{
		return set == null ? IntUnaryOperator.identity() : Indices.elementToIndexMap(set);
	}

	public static IntIterator iterator(MEC mec, int length)
	{
		return mec == null ? IntIterators.fromTo(0, length) : mec.states.iterator();
	}

	public static IntIterator iterator(IntCollection set, int length)
	{
		return set == null ? IntIterators.fromTo(0, length) : set.iterator();
	}
}
