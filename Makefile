init:
	pip install -r requirements.txt

test:
	pytest --cov=flametransfer tests

.PHONY: init test
