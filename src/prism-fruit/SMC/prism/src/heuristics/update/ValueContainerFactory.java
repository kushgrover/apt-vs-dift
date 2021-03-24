package heuristics.update;

import heuristics.search.HeuristicFactory;

public final class ValueContainerFactory
{
	private ValueContainerFactory() {}

	public static StateValueContainer createContainer(HeuristicFactory.HeuristicType heuristicType, boolean parallel) {
		final StateValueContainer container;
		if (heuristicType == HeuristicFactory.HeuristicType.RTDP_UNBOUNDED) {
			container = new DefaultStateValueContainer();
		} else {
			container = new DefaultDepthStateValueContainer();
		}
		if (parallel) {
			return new ThreadSafeStateValueContainer(container);
		}
		return container;
	}
}
