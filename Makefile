.PHONY: test formal release

test:
	./scripts/run_tests.sh

formal:
	./scripts/check_formal.sh

release:
	.venv/bin/python release/hadamard_7572_announcement/prepare.py
