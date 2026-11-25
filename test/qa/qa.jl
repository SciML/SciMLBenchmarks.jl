using JET
using SciMLBenchmarks
using SciMLTesting

run_qa(
    SciMLBenchmarks;
    # Preserving the sweep thread count requires the startup options; the public API gives only the total.
    ei_kwargs = (; all_qualified_accesses_are_public = (; ignore = (:JLOptions,)))
)
