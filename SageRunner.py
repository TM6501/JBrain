# This file is meant to be used from within C++ and is thus written
# to be as easy-to-use as possible from the clunky C++ Python API.
# As much as possible is hard-coded.
# import sys

# Because the PyWideStringList_Append function (C++) doesn't actually
# append, but rather replaces, we need to put the whole path back in place
# ourselves, here:
# sys.path = ['', 'C:\\Program Files\\Python310\\python310.zip',
#          'C:\\Program Files\\Python310\\DLLs',
#		   'C:\\Program Files\\Python310\\lib',
#		   'C:\\Program Files\\Python310',
#		   'C:\\Users\\josep\\AppData\\Roaming\\Python\\Python310\\site-packages',
#		   'C:\\Program Files\\Python310\\lib\\site-packages',
#		   'C:\\Users\\josep\\AppData\\Roaming\\Python\\Python310\\Scripts',
#		   'C:\\Users\\josep\\AppData\\Roaming\\Python\\Python310\\site-packages\\torch\\lib',
#		   'C:\\Users\\josep\\AppData\\Roaming\\Python\\Library\\bin']

import torch
import math
import numpy as np
import gymnasium as gym
from pathlib import Path

from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env, DummyVecEnv
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold

# Global variables to make life easy:
gModel = []
gEnv = None

# We want a simple 1-environment agent:
def getEnv(renderMode='None'):
	return DummyVecEnv([lambda: gym.make("CartPole-v1", render_mode=renderMode)])

# Training is independent and won't be called from C++:
def train(saveLocation):
	vec_env = getEnv()
	eval_env = gym.make("CartPole-v1")
	
	callback_on_best = StopTrainingOnRewardThreshold(reward_threshold=500, verbose=1)
	eval_callback = EvalCallback(eval_env, callback_on_new_best=callback_on_best, verbose=1)	

	model = PPO("MlpPolicy", vec_env, verbose=1)

	model.learn(total_timesteps=100000, callback=eval_callback)
	model.save(saveLocation)

def loadMultiplePPOSages(directory : str):
	global gModel
	searchDir = Path(directory)

	# Assume every zip file in the folder is a PPO model:
	for x in searchDir.iterdir():
	    if str(x).endswith(".zip"):
		    gModel.append(PPO.load(str(x)[:-4])) # Cut off the ".zip"

def load():
	global gModel
	global gEnv

	# If there is a folder of sages, load them. If not, load a single sage:
	if Path("sageCommittee").is_dir():
		loadMultiplePPOSages("sageCommittee")
	else:
		gModel = [PPO.load('ppo_sage')]

	gEnv = getEnv(renderMode='None')

def reset():
	global gEnv
	obs = gEnv.reset()
	return obs[0]  # Vector environment headaches.

# Reminder, as a sage it wants to return an observation and what the sage-model
# would do under those circumstances.
def step(actualAction: list):
    global gEnv
    global gModel

    inAction = np.array([actualAction])	
    # inAction = actualAction
    obs, rewards, dones, info = gEnv.step(inAction)
    actions = getModelAction(obs)

	# Vec-env returns need to be converted to what the singular
	# agent needs:
    return obs[0], rewards[0], dones[0], actions[0]

def getModelAction(obs):
    global gEnv
    global gModel

	# Get the first action:
    action, _ = gModel[0].predict(obs)

	# Add in the rest (possibly none):
    for i in range(1, len(gModel), 1):
        tmpAction, _ = gModel[i].predict(obs)
        action[0] += tmpAction[0]

	# Get the average of the actions as the vote outcome:
    action[0] = int(round(action[0] / len(gModel), 0))

    return action

def test():
    vec_env = getEnv(renderMode='human')
    models = load()

    obs = vec_env.reset()
    maxSteps = 20001
    # With a simple environment like this, what we want is in subscript
    # 0 of each of action, rewards, dones, and info.
    for i in range(1, maxSteps, 1):
        action = getModelAction(obs)
        obs, rewards, dones, info = vec_env.step(action)
        vec_env.render("human")
        if dones[0]:
            print(f"Step {i}: Done. Resetting the environment.")

if __name__ == "__main__":
    # for i in range(7):
        # print(f"Training model {i}...")
        # train(f"sageCommittee\\ppo_sage_{i}")
    test()
