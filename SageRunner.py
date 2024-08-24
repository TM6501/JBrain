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

from stable_baselines3 import PPO
from stable_baselines3.common.env_util import make_vec_env, DummyVecEnv
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold

# Global variables to make life easy:
gModel = None
gEnv = None

# We want a simple 1-environment agent:
def getEnv(renderMode='None'):
	return DummyVecEnv([lambda: gym.make("CartPole-v1", render_mode=renderMode)])

# Training is independent and won't be called from C++:
def train():
	vec_env = getEnv()
	eval_env = gym.make("CartPole-v1")
	
	callback_on_best = StopTrainingOnRewardThreshold(reward_threshold=500, verbose=1)
	eval_callback = EvalCallback(eval_env, callback_on_new_best=callback_on_best, verbose=1)	

	model = PPO("MlpPolicy", vec_env, verbose=1)

	model.learn(total_timesteps=100000, callback=eval_callback)
	model.save("ppo_sage")

def load():
	global gModel
	global gEnv

	gModel = PPO.load('ppo_sage')
	gEnv = genEnv(renderMode='None')

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
	obs, rewards, dones, info = gEnv.step(inAction)
	actions, states = gModel.predict(obs)

	# Vec-env returns need to be converted to what the singular
	# agent needs:
	return obs[0], rewards[0], dones[0], actions[0]

def test():
	vec_env = getEnv(renderMode='human')
	model = PPO.load("ppo_sage")

	obs = vec_env.reset()
	maxSteps = 2000
	# With a simple environment like this, what we want is in subscript
	# 0 of each of action, rewards, dones, and info.
	for i in range(maxSteps):
		action, _states = model.predict(obs)
		obs, rewards, dones, info = vec_env.step(action)
		vec_env.render("human")
		if dones[0]:
			print(f"Step {i}: Done. Resetting the environment.")


if __name__ == "__main__":
	# train()
	test()