ENV_DIR := /home/oneai/envs/

PY_ENV_NAME := enact_py_env

PY_ENV_PATH := $(ENV_DIR)$(PY_ENV_NAME)

CONFIG_PATH ?= config/configs.yaml

create-env:
	conda create --prefix $(PY_ENV_PATH) python=3.9

run_enact:
	bash setup_py_env.sh $(PY_ENV_PATH)
	bash run_enact.sh $(PY_ENV_PATH) ${CONFIG_PATH}

setup_py_env:
	bash setup_py_env.sh $(PY_ENV_PATH)

run_cell_ann_eval:
	bash setup_py_env.sh $(PY_ENV_PATH)
	bash run_cell_ann_eval.sh $(PY_ENV_PATH) 

reproduce_results:
	bash setup_py_env.sh $(PY_ENV_PATH)
	bash reproduce_paper_results.sh $(PY_ENV_PATH)